#' Create a lst of train/test indices
#'
#' @param n number of rows in data set.
#' @param nfold number of folds for cross-validation.
#' @importFrom purrr map
#' @importFrom caret createFolds
#'
#' @export
#' @keywords internal
get_cv_indices <- function(y, nfold = 4) {

  ## set up CV
  ind_test <- createFolds(y, k = nfold)
  ind_ret <- map(
    ind_test,
    ~ list(
      ind_train = setdiff(seq_along(y), .x),
      ind_test  = .x
    )
  )

}

#' Split data into training and eval data for watchlist
#'
#' @export
#' @keywords internal
get_split_indices <- function(n, split_frac = .7) {

  n_train   <- round(n  * split_frac, 0)
  train_idx <- sample(seq_len(n), n_train)
  eval_idx  <- setdiff(seq_len(n), train_idx)

  list(ind_train = train_idx, ind_eval = eval_idx)

}

#' Draw random parameters from prespecified parameter space
#'
#' Checks if integer or double. Samples integers or draws from uniform distribution
#' respectively.
#' @export
random_params <- function(param_space) {

  purrr::map(
    param_space,
    ~ {
        if(is.integer(.x[1])) {
          sample(seq(.x[1], .x[2], by = 1L), 1)
        } else {
          runif(1, .x[1], .x[2])
        }
      }
    )
}
