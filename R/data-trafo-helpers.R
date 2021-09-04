#' Transform data suitable for xgboost applications
#' @export
as_xgb_data <- function(x, ...) {

  UseMethod("as_xgb_data", x)

}


#' @rdname as_xgb_data
#' @importFrom stats model.matrix
#' @inherit as_xgb_data
#TODO: remove hard coded vars (e.g. offset)
as_xgb_data.data.frame <- function(x, features, label, set_margin = TRUE, ...) {

  mm <- model.matrix(~ ., x[, features, drop = FALSE])[, -1]
  xgb_data <- xgb.DMatrix(
    data = mm,
    label = as.numeric(x[[label]])
  )
  if(set_margin) {
    setinfo(xgb_data, "base_margin", x[["offset"]])
  }

  xgb_data

}

#' @rdname as_xgb_data
#' @inherit as_xgb_data
# TODO: remove hard coded vars (should be stored in ped object)
as_xgb_data.ped <- function(x, set_margin = TRUE, ...) {

  omit_vars <- attributes(x)[c("id_var", "intvars")] %>% unlist()
  keep_vars <- unique(c("tend", setdiff(colnames(x), omit_vars)))

  as_xgb_data(
    as.data.frame(x),
    features   = keep_vars,
    label      = "ped_status",
    set_margin = set_margin)

}

#' Transform new data to xgb-PED format using same arguments as previous transformation
#'
#' @keywords internal
as_xgb_newdata <- function(x, newdata, ...) {
  UseMethod("as_xgb_newdata", x)
}

#' @inherit as_xgb_newdata
#' @keywords internal
as_xgb_newdata.ped <- function(x, newdata, set_margin = TRUE, ...) {

    ped <- as_ped(x, newdata)
    vars <- setdiff(attr(x, "names"), attr(x, "intvars"))

    as_xgb_data(
      ped,
      features   = vars,
      label      = "ped_status",
      set_margin = set_margin)


}
#' @inherit as_xgb_newdata
#' @keywords internal
as_xgb_newdata.pam_xgb <- function(x, newdata, ...) {

    ped_newdata <- as_ped(x, newdata)

    as_xgb_data(ped_newdata, ...)

}

#' Transform newdata to PED format based on fitted xgb model
#'
#' @rdname as_ped
#' @importFrom pammtools as_ped split_data is.ped
#' @export
as_ped.pam_xgb <- function(data, newdata, ...) {

  if (is.ped(newdata)) {
    stop("newdata already in ped format.")
  }

  trafo_args      <- data[["trafo_args"]]
  trafo_args$data <- newdata
  do.call(pammtools::as_ped, trafo_args)

}

#' @rdname as_ped
as_ped.ped <- function(data, newdata, ...) {

  if (pammtools::is.ped(newdata)) {
    stop("newdata already in ped format.")
  }

  trafo_args <- attr(data, "trafo_args")
  trafo_args[["data"]] <- newdata
  do.call(pammtools::as_ped, trafo_args)

}
