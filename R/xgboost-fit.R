#' Fit xgboost model to piece-wise exponential data
#'
#' Transforms a data set of class \code{ped} to a format suitable for
#' \code{xgboost}, then runs the xgboost model. All parameters available in
#' \code{\link[xgboost]{xgboost}} can be passed on via ellipsis.
#' The \code{base_score} is set to \code{1}.
#' @importFrom xgboost xgboost xgb.DMatrix setinfo
#' @param data An objet of class \code{ped}.
#' @param ... Further arguments passed to \code{\link[xgboost]{xgboost}}.

#' @export
xgboost.ped <- function(data = NULL, ...) {

  attr_ped <- attributes(data)
  trafo_args <- attr_ped[["trafo_args"]]
  attr_ped[["trafo_args"]] <- NULL

  xgb_data <- as_xgb_data(data)

  pam_xgb <- xgboost(
    xgb_data,
    ...,
    base_score = 1,
    objective = "count:poisson")
  omit_vars <- attributes(data)[c("id_var", "intvars")] %>% unlist()
  keep_vars <- c("tend", setdiff(colnames(data), omit_vars))
  attr(pam_xgb, "attr_ped") <- attr_ped
  pam_xgb[["orig_features"]] <- keep_vars
  pam_xgb[["trafo_args"]] <- trafo_args
  class(pam_xgb) <- c("pam_xgb", class(pam_xgb))

  pam_xgb

}

#' PAMM wrapper for xgb.train
#'
#' @inheritParams xgboost::xgb.train
#' @importFrom xgboost xgb.train
#' @importFrom pammtools as_ped
#' @export
xgb.train.ped <- function(
  params                = list(),
  data,
  nrounds,
  watchlist             = list(),
  obj                   = NULL,
  feval                 = NULL,
  verbose               = 1,
  print_every_n         = 1L,
  early_stopping_rounds = NULL,
  maximize              = NULL,
  save_period           = NULL,
  save_name             = "xgboost.model",
  xgb_model             = NULL,
  callbacks             = list(),
  nthread               = 1L,
  base_score            = 1,
  ...
) {

  attr_ped                 <- attributes(data)
  trafo_args               <- attr_ped[["trafo_args"]]
  attr_ped[["trafo_args"]] <- attr_ped[["row_names"]] <-
    attr_ped[["class"]] <- NULL

  # xgb_data <- as_xgb_data()


  pam_xgb <- xgb.train(
    params                = params,
    data                  = as_xgb_data(data),
    nrounds               = nrounds,
    watchlist             = map(watchlist, ~ as_xgb_data(as_ped(data, newdata = .x))),
    obj                   = obj,
    feval                 = feval,
    verbose               = verbose,
    print_every_n         = print_every_n,
    early_stopping_rounds = early_stopping_rounds,
    maximize              = maximize,
    save_period           = save_period,
    save_name             = save_name,
    xgb_model             = xgb_model,
    callbacks             = callbacks,
    base_score            = base_score,
    objective             = "count:poisson",
    nthread               = nthread,
    ...
  )

  omit_vars <- c(attributes(data)[c("id_var", "intvars")] %>% unlist(), "ped_status") %>%
    unique()
  keep_vars <- unique(c("tend", setdiff(colnames(data), omit_vars)))
  attr(pam_xgb, "attr_ped")  <- attr_ped
  pam_xgb[["trafo_args"]]    <- trafo_args
  pam_xgb[["orig_features"]] <- keep_vars
  class(pam_xgb)             <- c("pam_xgb", class(pam_xgb))

  pam_xgb

}


#' Runs xgb.train.ped on cross-validation sets
#'
#'
#' @inheritParams xgb.train.ped
#' @param nfold Number of cross-valdation folds.
#' @param ped_params List of parameters used to transform data into PED format.
#' @export
xgb.cv_pam <- function(
  params = list(),
  data,
  nrounds,
  nfold                 = 4,
  cv_indices,
  ped_params            = list(),
  nthread               = 1L,
  verbose               = FALSE,
  print_every_n         = 1L,
  early_stopping_rounds = NULL,
  ...
) {

  if (missing(cv_indices)) {
    cv_indices <- get_cv_indices(data$status, nfold)
  }

  cv_res <- map(
    cv_indices,
    ~{
      if (class(data) == "list") {
        data_train <- data[[1]][.x$ind_train, , drop = FALSE]
        data_test <- data[[1]][.x$ind_test, , drop =FALSE]
        ped_params[["data"]] <- list(data_train, data[[2]])
        watch <- list(data_test, data[[2]])
      } else {
        ped_params[["data"]] <- data[.x$ind_train, , drop = FALSE]
        watch <- data[.x$ind_test, , drop = FALSE]
      }
      # ped <- as_ped(data = ped_params$data, formula = ped_params$formula)
      ped <- do.call(as_ped, ped_params)
      xgb <- xgb.train.ped(
        params = params,
        data   = ped,
        watchlist = list(eval = watch),
        nthread = nthread,
        nrounds = nrounds,
        verbose = verbose,
        print_every_n = print_every_n,
        early_stopping_rounds = early_stopping_rounds)
    }
  )


}


#' Tune xgb pam
#'
#' @inheritParams xgb.cv_pam
#' @param param_df A data frame of parameter combinations to tune. One
#' row contains one parameter set that will be passed on to
#' \code{params} in \code{xgb.cv_pam}.
#' @importFrom prodlim Hist
#' @importFrom survival Surv
#' @importFrom pec pec crps
#' @export
xgb.tune_pam <- function(
  data,
  param_df,
  nrounds,
  cv_indices,
  nfold      = 4,
  ped_params = list(),
  nthread    = 1L,
  early_stopping_rounds = NULL,
  verbose = FALSE,
  print_every_n = 1L,
  ...
) {

  param_list <- split(param_df, seq_len(nrow(param_df)))

  if (missing(cv_indices)) {
    cv_indices <- get_cv_indices(data$status, nfold)
  }

  tune_res <- map_dfr(
    param_list,
    ~ {
      print(.x)
      cv_res <- xgb.cv_pam(
        data                  = data,
        params                = .x,
        nrounds               = nrounds,
        nfold                 = nfold,
        cv_indices            = cv_indices,
        ped_params            = ped_params,
        nthread               = nthread,
        early_stopping_rounds = early_stopping_rounds,
        verbose               = verbose,
        print_every_n         = print_every_n
      )
      cv_iter_best <- purrr::map_dbl(cv_res, ~.x$best_iteration)
      ibs_df <- purrr::map2_dfr(
        .x = cv_res,
        .y = cv_indices,
        ~ {
          test_dat <- data
          if(class(test_dat) == "list") {
            test_dat[[1]] <- data[[1]][.y$ind_test, , drop = FALSE]
          } else {
            test_dat <- test_dat[.y$ind_test, , drop = FALSE]
          }
          get_ibs(
            object  = .x,
            data    = test_dat,
            pec_params = list(
              formula = Hist(time, status) ~ 1,
              start   = .01,
              exact   = FALSE),
            keep_reference = FALSE)
        },
        .id = "fold"
      )

      ibs_df %>%
        group_by(fold, cause) %>%
        summarize(IBS = mean(IBS)) %>%
        group_by(fold) %>%
        summarize(IBS = mean(IBS)) %>%
        mutate(nrounds = cv_iter_best)

    },
    .id = "id"
  )

  tune_res %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(
      IBS = mean(IBS),
      nrounds = round(mean(nrounds))) %>%
    dplyr::mutate(param_set = param_list)

}


#' Runs xgb.train.ped on cross-validation sets
#'
#'
#' @inheritParams xgb.train.ped
#' @param nfold Number of cross-valdation folds.
#' @param ped_params List of parameters used to transform data into PED format.
#' @export
rf.cv_pam <- function(
  params = list(),
  data,
  nrounds = 1L,
  nfold                 = 4,
  cv_indices,
  ped_params            = list(),
  nthread               = 1L,
  verbose               = FALSE,
  print_every_n         = 1L,
  early_stopping_rounds = NULL,
  ...
) {

  if (missing(cv_indices)) {
    cv_indices <- get_cv_indices(data$status, nfold)
  }

  cv_res <- map(
    cv_indices,
    ~{
      if (class(data) == "list") {
        data_train <- data[[1]]
        ped_params[["data"]] <- list(data_train, data[[2]])
      } else {
        ped_params[["data"]] <- data
      }
      # ped <- as_ped(data = ped_params$data, formula = ped_params$formula)
      ped <- do.call(as_ped, ped_params)
      xgb <- xgb.train.ped(
        params = params,
        data   = ped,
        nthread = nthread,
        nrounds = nrounds,
        verbose = verbose,
        print_every_n = print_every_n)
    }
  )


}


#' Tune xgb pam
#'
#' @inheritParams rf.cv_pam
#' @param param_df A data frame of parameter combinations to tune. One
#' row contains one parameter set that will be passed on to
#' \code{params} in \code{xgb.cv_pam}.
#' @importFrom prodlim Hist
#' @importFrom survival Surv
#' @importFrom pec pec crps
#' @export
rf.tune_pam <- function(
  data,
  param_df,
  nrounds = 1L,
  cv_indices,
  nfold      = 4,
  ped_params = list(),
  nthread    = 1L,
  early_stopping_rounds = NULL,
  verbose = FALSE,
  print_every_n = 1L,
  ...
) {

  param_list <- split(param_df, seq_len(nrow(param_df)))

  if (missing(cv_indices)) {
    cv_indices <- get_cv_indices(data$status, nfold)
  }

  tune_res <- map_dfr(
    param_list,
    ~ {
      print(.x)
      cv_res <- xgb.cv_pam(
        data                  = data,
        params                = .x,
        nrounds               = nrounds,
        nfold                 = nfold,
        cv_indices            = cv_indices,
        ped_params            = ped_params,
        nthread               = nthread,
        early_stopping_rounds = early_stopping_rounds,
        verbose               = verbose,
        print_every_n         = print_every_n
      )
      ibs_df <- purrr::map2_dfr(
        .x = cv_res,
        .y = cv_indices,
        ~ {
          test_dat <- data
          if(class(test_dat) == "list") {
            test_dat[[1]] <- data[[1]][.y$ind_test, , drop = FALSE]
          } else {
            test_dat <- test_dat[.y$ind_test, , drop = FALSE]
          }
          get_ibs(
            object  = .x,
            data    = test_dat,
            pec_params = list(
              formula = Hist(time, status) ~ 1,
              start   = .01,
              exact   = FALSE),
            keep_reference = FALSE)
        },
        .id = "fold"
      )

      ibs_df %>%
        group_by(fold, cause) %>%
        summarize(IBS = mean(IBS)) %>%
        group_by(fold) %>%
        summarize(IBS = mean(IBS))

    },
    .id = "id"
  )

  tune_res %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(
      IBS = mean(IBS),
      nrounds = round(mean(nrounds))) %>%
    dplyr::mutate(param_set = param_list)

}
