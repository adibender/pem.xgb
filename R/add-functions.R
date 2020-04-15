#' Add predicted (cumulative) hazard to data set
#'
#' Add (cumulative) hazard based on the provided data set and model.
#'
#'
#' @rdname add_hazard2
#' @inheritParams pammtools::add_hazard
#' @inheritParams xgboost::xgboost
#' @param ... Further arguments passed to \code{\link[xgboost]{predict.xgb.Booster}}.
#' @param overwrite Should hazard columns be overwritten if already present in
#' the data set? Defaults to \code{FALSE}.
#' @import checkmate dplyr
#' @importFrom stats predict
#' @export
add_hazard2 <- function(
  newdata,
  object,
  reference = NULL,
  overwrite = FALSE,
  time_var  = NULL,
  ...)  {

  if (!overwrite) {
    if ("hazard" %in% names(newdata)) {
      stop("Data set already contains 'hazard' column.
        Set `overwrite=TRUE` to overwrite")
    }
  } else {
      rm.vars <- intersect(
        c("hazard", "se", "ci_lower", "ci_upper"),
        names(newdata))
      newdata <- newdata %>% select(-one_of(rm.vars))
  }


  get_hazard2(object, newdata, reference = reference, ...)

}


get_hazard2 <- function(object, newdata, reference = NULL, ...) {
  UseMethod("get_hazard2", object)
}

get_hazard2.xgb.Booster <- function(
  object,
  newdata,
  reference = NULL,
  time_var = NULL,
  ...) {

  assert_data_frame(newdata, all.missing = FALSE)

  if (is.null(time_var)) {
    time_var <- "tend"
  } else {
    assert_string(time_var)
    assert_choice(time_var, colnames(newdata))
  }

  features <- unique(object$orig_features)
  X <- as_xgb_data(newdata, features = features, label = "ped_status",
    set_margin = FALSE)
  hazard_name <- paste0("hazard_", class(object)[1])
  newdata[[hazard_name]] <- xgboost:::predict.xgb.Booster(object, X)

  newdata

}

preproc_reference <- function(reference, cnames, n_rows) {

  # check that provided variables contained in newdata
  names_ref <- names(reference)
  if (!check_subset(names_ref, cnames)) {
    stop(paste0("Columns in 'reference' but not in 'newdata':",
      paste0(setdiff(names_ref, cnames), collapse = ",")))
  }
  # transform to list if inherits from data frame, so it can be processed
  # in mutate via !!!
  if (inherits(reference, "data.frame")) {
    if (!(nrow(reference) == n_rows | nrow(reference) == 1)) {
      stop("If reference is provided as data frame, number of rows must be
        either 1 or the number of rows in newdata.")
    }
    reference <- as.list(reference)
  }

  reference

}



#' @rdname add_hazard2
#' @inheritParams add_hazard2
#' @param interval_length The variable in newdata containing the interval lengths.
#' Can be either bare unquoted variable name or character. Defaults to \code{"intlen"}.
#' @importFrom dplyr bind_cols
#' @export
add_cumu_hazard2 <- function(
  newdata,
  object,
  overwrite = FALSE,
  time_var  = NULL,
  interval_length = "intlen",
  ...)  {

  interval_length <- quo_name(enquo(interval_length))

  if (!overwrite) {
    if ("cumu_hazard" %in% names(newdata)) {
      stop(
        "Data set already contains 'hazard' column.
        Set `overwrite=TRUE` to overwrite")
    }
  } else {
      rm.vars <- intersect(c("cumu_hazard", "cumu_lower", "cumu_upper"),
        names(newdata))
      newdata <- newdata %>% select(-one_of(rm.vars))
  }

  get_cumu_hazard2(object, newdata, time_var = time_var,
    interval_length = interval_length, ...)

}

#' Calculate cumulative hazard
#'
#' @inheritParams add_cumu_hazard2
#' @import checkmate dplyr
#' @importFrom rlang UQ sym quo_name .data .env
#' @importFrom purrr map_lgl
#' @keywords internal
get_cumu_hazard2 <- function(
  object,
  newdata,
  time_var   = NULL,
  interval_length = "intlen",
  ...)  {

  assert_character(interval_length)
  assert_subset(interval_length, colnames(newdata))
  assert_data_frame(newdata, all.missing = FALSE)

  interval_length <- sym(interval_length)
  hazard_name <- paste0("hazard_", class(object)[1])
  cumu_name <- paste0("cumu_hazard_", class(object)[1])

  mutate_args  <- list(cumu_hazard = quo(cumsum(.data[[hazard_name]] *
    (!!interval_length))))
  haz_vars_in_data <- map(c(hazard_name, "se"),
    ~ grep(.x, colnames(newdata), value = TRUE, fixed = TRUE)) %>% flatten_chr()
  vars_exclude <- c(hazard_name)

  newdata <- get_hazard2(object, newdata, time_var = time_var, ...)

  newdata <- newdata %>%
    mutate(!!!mutate_args)
  newdata <- newdata %>% rename({{cumu_name}} := "cumu_hazard")

  vars_exclude <- setdiff(vars_exclude, haz_vars_in_data)
  if (length(vars_exclude) != 0 ) {
    newdata <- newdata %>% select(-one_of(vars_exclude))
  }

  newdata

}


#' Add survival probability estimates
#'
#' Given suitable data (i.e. data with all columns used for estimation of the model),
#' this functions adds a column \code{surv_prob} containing survival probabilities
#' for the specified covariate and follow-up information.
#'
#' @inherit add_cumu_hazard2
#' @export
add_surv_prob2 <- function(
  newdata,
  object,
  overwrite  = FALSE,
  time_var   = NULL,
  interval_length = "intlen",
  ...)  {

  interval_length <- quo_name(enquo(interval_length))

  if (!overwrite) {
    if ("surv_prob" %in% names(newdata)) {
      stop("Data set already contains 'surv_prob' column.
        Set `overwrite=TRUE` to overwrite")
    }
  } else {
      rm.vars <- intersect(c("surv_prob"), names(newdata))
      newdata <- newdata %>% select(-one_of(rm.vars))
  }

  get_surv_prob2(object, newdata,
    time_var = time_var, interval_length = interval_length, ...)

}


#' Calculate survival probabilities
#'
#' @inheritParams add_surv_prob2
#' @importFrom purrr flatten_chr
#' @keywords internal
get_surv_prob2 <- function(
  object,
  newdata,
  time_var        = NULL,
  interval_length = "intlen",
  ...) {

  assert_character(interval_length)
  assert_subset(interval_length, colnames(newdata))
  assert_data_frame(newdata, all.missing = FALSE)

  interval_length <- sym(interval_length)
  hazard_name <- paste0("hazard_", class(object)[1])
  surv_name <- paste0("surv_prob_", class(object)[1])

  mutate_args  <- list(surv_prob = quo(exp(-cumsum(.data[[hazard_name]] *
    (!!interval_length)))))
  haz_vars_in_data <- map(c(hazard_name),
    ~grep(.x, colnames(newdata), value = TRUE, fixed = TRUE)) %>% flatten_chr()
  vars_exclude <- c(hazard_name)

  newdata <- get_hazard2(
      object = object,
      newdata,
      time_var = time_var)

  newdata <- newdata %>%
    mutate(!!!mutate_args)
  newdata <- newdata %>% rename({{surv_name}} := "surv_prob")

  vars_exclude <- setdiff(vars_exclude, haz_vars_in_data)
  if (length(vars_exclude) != 0 ) {
    newdata <- newdata %>% select(-one_of(vars_exclude))
  }

  newdata

}


# add_cif <- function() {
#   newdata,
#   object,
#   overwrite  = FALSE,
#   time_var   = NULL,
#   interval_length = "intlen",
#   ...)  {

#   interval_length <- quo_name(enquo(interval_length))

#   if (!overwrite) {
#     if ("surv_prob" %in% names(newdata)) {
#       stop("Data set already contains 'surv_prob' column.
#         Set `overwrite=TRUE` to overwrite")
#     }
#   } else {
#       rm.vars <- intersect(c("surv_prob"), names(newdata))
#       newdata <- newdata %>% select(-one_of(rm.vars))
#   }

#   get_surv_prob2(object, newdata,
#     time_var = time_var, interval_length = interval_length, ...)

# }

# #' Calculate all-cause survival probabilities
# #'
# #' @inheritParams add_surv_prob2
# #' @importFrom purrr flatten_chr
# #' @keywords internal
# get_all_cause_surv_prob <- function(
#   object,
#   newdata,
#   time_var        = NULL,
#   interval_length = "intlen",
#   status_var      = "type",
#   ...) {

#   assert_character(interval_length)
#   assert_subset(interval_length, colnames(newdata))
#   assert_data_frame(newdata, all.missing = FALSE)


#   interval_length <- sym(interval_length)
#   hazard_name <- paste0("hazard_", class(object)[1])
#   hazard_names <- paste0(hazard_name, seq_len(unique(newdata[[status_var]])))
#   surv_name <- paste0("surv_prob_allcause", class(object)[1])

#   mutate_args  <- list(surv_prob = quo(exp(-cumsum(.data[[hazard_name]] *
#     (!!interval_length)))))
#   haz_vars_in_data <- map(c(hazard_name),
#     ~grep(.x, colnames(newdata), value = TRUE, fixed = TRUE)) %>% flatten_chr()
#   vars_exclude <- c(hazard_name)

#   newdata1 <- get_hazard2(
#       object = object,
#       dplyr::mutate(newdata$status = status_vals[1],
#       time_var = time_var)

#   newdata <- newdata %>%
#     mutate(!!!mutate_args)
#   newdata <- newdata %>% rename({{surv_name}} := "surv_prob")

#   vars_exclude <- setdiff(vars_exclude, haz_vars_in_data)
#   if (length(vars_exclude) != 0 ) {
#     newdata <- newdata %>% select(-one_of(vars_exclude))
#   }

#   newdata

# }
