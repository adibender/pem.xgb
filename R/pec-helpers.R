#' Extract IBS in tidy format
#'
#' @inheritParams pec::pec
#' @inheritParams pec::crps
#' @importFrom stats quantile
#' @param keep_reference Logical. Should Brier Score of Reference model be kept
#' in results table.
#' @export
get_ibs <- function(
  object,
  data,
  pec_params = list(
    formula = Surv(time, status) ~ 1,
    exact   = FALSE,
    start   = .01),
  keep_reference = TRUE,
  q_last = .8,
  q_eval = c(.25, .5, .75),
  ...) {

  if (class(data)[1] == "list") {
    status_values <- sort(unique(data[[1]]$status))
  } else {
    status_values <- sort(unique(data$status))
  }

  cause_values <- status_values[status_values != 0]

  res <- map_dfr(
    cause_values,
    ~ {
      times <- times_list(data, q_last = q_last, q_eval = q_eval,
          status_value = .x)
      time_seq <- seq(.01, times$t_last, length.out = 500L)
      if (class(data)[1] == "list") {
        if(class(object)[1] == "list") {
          pred <- predictSurvProb(object[[1]], data, times = time_seq)
        } else {
          pred <- predictSurvProb(object, data, times = time_seq)
        }
      }
      pec_params$object <- if (class(data)[1] == "list") {list("pam_xgb" = pred)} else {object}
      pec_params$data <- if (class(data)[1] == "list") { data[[1]]} else{ data }
      pec_params$times <- time_seq
      pec_params$cause <- .x

      if(length(cause_values) > 1) {
        pec_params$formula <- update(pec_params$formula, Hist(time, status)~.)
      }

      pec <- do.call(pec::pec, pec_params)

      ibs <- pec::crps(pec, times = times$eval_times) %>%
        as.data.frame() %>%
        group_by(method) %>%
        mutate(
          time     = times$eval_times,
          quantile = q_eval) %>%
        ungroup() %>%
        mutate(cause = .x)

      ibs

    }
  )

  if (!keep_reference) {
    res <- res %>%
      dplyr::filter(method != "Reference")
  }

  res

}


#' Look up Quantiles of event times in data
#'
#' @keywords internal
times_list <- function(data, q_last = .8, q_eval = c(.25, .5, .75),
  time_var = "time", status_var = "status", status_value = 1) {

  if (class(data)[1] == "list") {
    data <- data[[1]]
  }

  t_last     <- quantile(data[[time_var]][data[[status_var]] == status_value], q_last)
  eval_times <- quantile(data[[time_var]][data[[status_var]] == status_value], q_eval)

  list(t_last = t_last, eval_times = eval_times)

}


#' C-Index from pec package
#'
get_cindex <- function(
  object,
  data,
  cind_params = list(
    formula = Surv(time, status) ~ 1,
    exact   = FALSE,
    start   = .01),
  eval_times = NULL,
  q_last = .8,
  q_eval = c(.25, .5, .75)) {

  if(is.null(eval_times)) {
    times <- times_list(data, q_last = q_last, q_eval = q_eval)
    eval_times <- times$eval_times
  }


  cind_params$object <- object
  cind_params$data <- data
  cind_params$eval.times <- eval_times
  cind <- do.call(pec::cindex, cind_params)

  cind$AppCindex %>% as.data.frame() %>%
    mutate(eval_time = eval_times) %>%
    tidyr::gather("method", "cindex", -eval_time)

}


#' Transform crps object to data.frame
#'
#' A\code{as.data.frame} S3 method for objects of class \code{\link[pec]{crps}}.
#'
#' @inheritParams base::as.data.frame
#' @param x An object of class \code{crps}. See \code{\link[pec]{crps}}.
#' @importFrom tidyr pivot_longer
#'
#' @export
as.data.frame.crps <- function(x, row.names = NULL, optional = FALSE, ...) {

  m <- matrix(x, nrow = dim(x)[1], ncol = dim(x)[2])
  colnames(m) <- attr(x, "dimnames")[[2]]

  m <- as.data.frame(m)
  m$method <- attr(x, "dimnames")[[1]]

  m <- m %>%
    pivot_longer(cols = -.data$method, values_to = "IBS") %>%
    dplyr::rename(time = .data$name)

}
