context("Test fitting xgb models for survival data")


test_that("xgb models are fit to PED data", {

  data("tumor", package = "pammtools")
  ped <- pammtools::as_ped(tumor[1:100,], Surv(days, status) ~ . )

  ## xgboost
  mod <- xgboost.ped(ped, nrounds = 100, print_every_n = 50)
  expect_class(mod, c("pam_xgb", "xgb.Booster"))
  expect_class(mod[["trafo_args"]], "list")
  expect_identical(length(mod[["trafo_args"]]), 3L)
  expect_identical(names(mod[["trafo_args"]]), c("formula", "cut", "id"))

  ## xgb.train
  mod2 <- xgb.train.ped(
    params    = list(
      max_depth = c(3, 5),
      eta       = 0.3
    ),
    data = as_ped(tumor[1:100,], Surv(days, status)~.),
    nrounds = 500L,
    watchlist = list(eval = tumor[201:300, ]
    ),
    verbose = FALSE,
    early_stopping_rounds = 50
  )

  expect_class(mod2, c("pam_xgb", "xgb.Booster"))
  expect_class(mod2[["trafo_args"]], "list")
  expect_identical(length(mod2[["trafo_args"]]), 3L)
  expect_identical(names(mod2[["trafo_args"]]), c("formula", "cut", "id"))


})
