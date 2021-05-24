#' Predict hazard, cumulative hazard or survival probability
#'
#' @param object An object of class \code{pam_xgb}
#' @param newdata A data set containing the same covariates as used for model
#' fitting. If of class \code{data.frame}, the function will try to transform
#' to the \code{xgb.DMatrix} format.
#' @param times A vector of times for which predictions should be generated
#' @param type The type of prediction desired. Either hazard (\code{type = "hazard"}),
#' cumulative hazard (\code{type = "cumu_hazard"}) or survival probability
#' (\code{type = "surv_prob"}).
#' @importFrom stats predict
#' @importFrom purrr map_dfr
#' @return A matrix of predictions containing one row per
#' observation (row in newdata) and 1 column per  specified time in the
#' \code{times} argument.
#' @export
predict.pam_xgb <- function(
  object,
  newdata,
  type = c("hazard", "cumu_hazard", "surv_prob"),
  set_margin = FALSE,
  ...) {

  type <- match.arg(type)
  brks <- attr(object, "attr_ped")$breaks
  # attr(object, "attr_ped") <- c(attr(object, "attr_ped"), status_var = "status")

  ped_newdata <- as_ped(object, newdata)
  vars <- setdiff(
    attr(ped_newdata, "names"),
    attr(ped_newdata, "intvars"))

  xgb_newdata   <- as_xgb_data(ped_newdata, set_margin = set_margin)

  ped_newdata[["pred"]] <- xgboost:::predict.xgb.Booster(
    object,
    xgb_newdata,
    ...)
  if (type == "cumu_hazard") {
    ped_newdata <- ped_newdata %>%
      group_by(.data$id) %>%
      mutate(pred = cumsum(.data$pred * exp(.data$offset)))
  }
  if (type == "surv_prob") {
     ped_newdata <- ped_newdata %>%
      group_by(.data[["id"]]) %>%
      mutate(pred = exp(-cumsum(.data$pred * exp(.data$offset))))
  }

  ped_newdata %>%
    group_by(.data[["id"]]) %>%
    filter(row_number() == n()) %>%
    pull(.data[["pred"]])

}

# check for time-dependent covariates.
has_tdc <- function(model) {
  any(c("ccr", "func") %in% names(attributes(model)[["attr_ped"]]))
}

get_new_ped <- function(object, newdata, times, attr_ped) {
  covars <- setdiff(attr_ped[["names"]], attr_ped[["intvars"]])
  if ("tend" %in% object$feature_names) {
    vars <- c("tend", covars)
  }
  id_var <- attr_ped[["id_var"]]
  brks <- attr_ped[["breaks"]]
  ped_times <- sort(unique(union(c(0, brks), times)))
  # extract relevant intervals only, keeps data small
  ped_times <- ped_times[ped_times <= max(times)]
  # obtain interval information
  ped_info <- get_intervals(brks, ped_times[-1])
  # add adjusted offset such that cumulative hazard and survival probability
  # can be calculated correctly
  ped_info[["offset"]] <- c(ped_info[["times"]][1], diff(ped_info[["times"]]))

  # create data set with interval/time + covariate info
  newdata[[id_var]] <- seq_len(nrow(newdata))
  new_ped <- pammtools::combine_df(ped_info, newdata[, c(id_var, covars)])
  new_ped$ped_status <- 1
  new_ped
}

get_all_intervals <- function(
  x,
  times,
  left.open        = TRUE,
  rightmost.closed = TRUE,
  ...) {

  # check inputs
  assert_numeric(times, lower = 0, finite = TRUE, all.missing = FALSE)

  int_df <- pammtools::int_info(x)
  int    <- findInterval(
    x                = times,
    vec              = sort(union(int_df$tstart, int_df$tend)),
    left.open        = left.open,
    rightmost.closed = rightmost.closed)

  #int_df %>%
  #  slice(int) %>%
  # Q: slice(int) will drop duplicated entries -- this is not what we want here, I think?
  # should this be changed in pammtools::get_intervals.default as well?
  int_df[int,] %>%
    mutate(times = times) %>%
    # arrange(times) %>%
    select(times, everything())

}

#' @importFrom tidyr fill
get_new_ped_tdc <- function(object, newdata, times, attr_ped) {
  # extract vars used in model fit
  covars <- setdiff(attr_ped[["names"]], attr_ped[["intvars"]])
  id_var <- attr_ped[["id_var"]]
  brks <- attr_ped[["breaks"]]
  # stopifnot(all(times <= max(brks)))

  #get time info for predictive TDCs
  ccr_tz_vars <- purrr::map_chr(attr_ped[["ccr"]][["ccr_list"]],
    ~.x[["tz_var"]]) %>% unique()
  # avoid assumption that newdata[[2]] is a single data.frame with TDCs on
  # the same time grid. could also be a list of data.frames (?)
  if (!is.data.frame(newdata[[2]])) {
    stop("multiple time scales for TDCs not implemented yet.")
  }
  #drop ids for which no time constant info is available:
  newdata[[2]] <- filter(newdata[[2]],
    !!sym(id_var) %in% newdata[[1]][[id_var]])

  tdc_times <- select(newdata[[2]], any_of(ccr_tz_vars)) %>%
    pull(1)

  ped_times <- sort(unique(c(0, brks, tdc_times, times)))
  # extract relevant intervals only, keeps data small
  ped_times <- ped_times[ped_times <= max(times)]
  # obtain interval information
  ped_info <- get_all_intervals(brks, ped_times[-1])
  # add adjusted offset such that cumulative hazard and survival probability
  # can be calculated correctly
  ped_info[["offset"]] <- c(ped_info[["times"]][1], diff(ped_info[["times"]]))

  # create data set with interval/time + time constant covariate info
  ccr_covars <- map(attr_ped[["ccr"]][["ccr_list"]],  ~.x[["col_vars"]]) %>%
    unlist() %>%  unique()
  const_covars <- setdiff(covars, ccr_covars)
  new_ped_const <- pammtools::combine_df(ped_info,
    newdata[[1]][, c(id_var, const_covars)])
  new_ped_const$ped_status <- 1

  # create data set with time dependent covariates (last value carried forward)
  suppressMessages({
  new_ped_tdc <- rename(newdata[[2]], times = sym(ccr_tz_vars)) %>%
    full_join(select(new_ped_const, all_of(c(id_var, "times")))) %>%
    filter(.data$times <= max(ped_times)) %>%
    arrange(!!sym(id_var), !!sym("times")) %>%
    tidyr::fill(all_of(ccr_covars)) %>%
    filter(.data$times != min(ped_times)) #drop rows for "time origin"

  inner_join(new_ped_const, new_ped_tdc)

  })

}
#' S3 method for compatibility with package pec
#'
#' @importFrom pec predictSurvProb
#' @importFrom purrr map
#' @importFrom pammtools get_intervals
#' @export
predictSurvProb.pam_xgb <- function(
    object,
    newdata,
    times) {

  attr_ped <- attributes(object)[["attr_ped"]]
  covars <- setdiff(attr_ped[["names"]], attr_ped[["intvars"]])
  if ("tend" %in% object$feature_names) {
    vars <- c("tend", covars)
  }
  id_var <- attr_ped[["id_var"]]

  if (has_tdc(object)) {
    new_ped <- get_new_ped_tdc(object, newdata, times, attr_ped)
  } else {
    new_ped <- get_new_ped(object, newdata, times, attr_ped)
  }

  new_ped[["pred"]] <- xgboost:::predict.xgb.Booster(
    object,
    as_xgb_data(new_ped, vars, label = "ped_status", set_margin = FALSE))## Note label not needed, set any
  new_ped <- new_ped %>%
    arrange(.data$id, .data$times) %>%
    group_by(.data$id) %>%
    mutate(pred = exp(-cumsum(.data$pred * .data$offset))) %>%
    ungroup() %>%
    filter(.data[["times"]] %in% .env[["times"]])

  id <- unique(new_ped[[id_var]])
  pred_list <- map(
    id,
    ~ new_ped[new_ped[[id_var]] == .x, "pred"] %>% pull("pred"))

  do.call(rbind, pred_list)

}


#' S3 method for compatibility with package pec
#'
#' This function is needed to use \code{pec::pec} in the competing risks setting.
#'
#' @importFrom pec predictEventProb
#' @importFrom purrr map
#' @importFrom pammtools get_intervals
#' @export
#' @rdname predictEventProb
predictEventProb.pam_xgb <- function(
  object,
  newdata,
  times,
  cause,
  ...
) {

  attr_ped <- attributes(object)[["attr_ped"]]
  id_var <- attr_ped[["id_var"]]
  brks <- attr_ped[["breaks"]]
  ped_times <- sort(unique(union(c(0, brks), times)))
  ped_times <- ped_times[ped_times <= max(times)]
  ped_info <- get_intervals(brks, ped_times[-1])
  ped_info[["offset"]] <- c(ped_info[["times"]][1], diff(ped_info[["times"]]))
  covars <- setdiff(attr_ped[["names"]], attr_ped[["intvars"]])
  if ("tend" %in% object$feature_names) {
    vars <- c("tend", covars)
  }
  newdata[["type"]] <- 1
  newdata[[id_var]] <- seq_len(nrow(newdata))
  new_ped <- pammtools::combine_df(ped_info, newdata[, c(id_var, covars)])
  new_ped$ped_status <- 1 # irrelevant
  new_ped$type <- 1
  new_ped[["csh1"]] <- xgboost:::predict.xgb.Booster(object,
        as_xgb_data(new_ped, vars, label = "ped_status", set_margin = FALSE))
  new_ped$type <- 2 # predict hazard for second cause
  new_ped[["csh2"]] <- xgboost:::predict.xgb.Booster(object,
        as_xgb_data(new_ped, vars, label = "ped_status", set_margin = FALSE))
  new_ped <- new_ped %>%
    arrange(.data$id, .data$times) %>%
    group_by(.data$id) %>%
    mutate(sp_all_cause = exp(-(
      cumsum(.data$csh1 * .data$offset) +
      cumsum(.data$csh2 * .data$offset)))) %>%
    mutate(cif1 = cumsum(.data$csh1 * (.data$sp_all_cause - 1e-5) * .data$offset)) %>%
    mutate(cif2 = cumsum(.data$csh2 * (.data$sp_all_cause - 1e-5) * .data$offset)) %>%
    ungroup() %>%
    filter(.data[["times"]] %in% .env[["times"]])
    id <- unique(new_ped[[id_var]])
    pred_list <- map(id, ~new_ped[new_ped[[id_var]] == .x,
        paste0("cif", cause)] %>%
        pull(paste0("cif", cause)))
    do.call(rbind, pred_list)

}


#' S3 method for compatibility with package pec
#'
#' This function is needed to use \code{pec::pec} in the competing risks setting.
#'
#' @importFrom pec predictEventProb
#' @importFrom purrr map
#' @importFrom pammtools get_intervals
#' @export
#' @rdname predictEventProb
predictEventProb.list <- function(
  object,
  newdata,
  times,
  cause,
  ...
) {

  attr_ped <- attributes(object[[1]])[["attr_ped"]]
  id_var <- attr_ped[["id_var"]]
  brks <- attr_ped[["breaks"]]
  ped_times <- sort(unique(union(c(0, brks), times)))
  ped_times <- ped_times[ped_times <= max(times)]
  ped_info <- get_intervals(brks, ped_times[-1])
  ped_info[["offset"]] <- c(ped_info[["times"]][1], diff(ped_info[["times"]]))
  covars <- setdiff(attr_ped[["names"]], attr_ped[["intvars"]])
  if ("tend" %in% object[[1]]$feature_names) {
    vars <- c("tend", covars)
  } else {
    vars <- covars
  }
  newdata[["type"]] <- 1
  newdata[[id_var]] <- seq_len(nrow(newdata))
  new_ped <- pammtools::combine_df(ped_info, newdata[, c(id_var, covars)])
  new_ped$ped_status <- 1 # irrelevant
  new_ped$type <- 1
  new_ped[["csh1"]] <- xgboost:::predict.xgb.Booster(object[[1]],
        as_xgb_data.data.frame(new_ped, vars, label = "ped_status", set_margin = FALSE))
  new_ped$type <- 2 # predict hazard for second cause
  new_ped[["csh2"]] <- xgboost:::predict.xgb.Booster(object[[2]],
        as_xgb_data.data.frame(new_ped, vars, label = "ped_status", set_margin = FALSE))
  new_ped <- new_ped %>%
    arrange(.data$id, .data$times) %>%
    group_by(.data$id) %>%
    mutate(sp_all_cause = exp(-(
      cumsum(.data$csh1 * .data$offset) +
      cumsum(.data$csh2 * .data$offset)))) %>%
    mutate(cif1 = cumsum(.data$csh1 * (.data$sp_all_cause - 1e-5) * .data$offset)) %>%
    mutate(cif2 = cumsum(.data$csh2 * (.data$sp_all_cause - 1e-5) * .data$offset)) %>%
    ungroup() %>%
    filter(.data[["times"]] %in% .env[["times"]])
    id <- unique(new_ped[[id_var]])
    pred_list <- map(id, ~new_ped[new_ped[[id_var]] == .x,
        paste0("cif", cause)] %>%
        pull(paste0("cif", cause)))
    do.call(rbind, pred_list)

}
