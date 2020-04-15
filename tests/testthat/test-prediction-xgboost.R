context("Test prediction functions for xgboost objects")

test_that("Generic predict function works", {
  library(pammtools)
  data("tumor", package = "pammtools")
  ped <- pammtools::as_ped(tumor[1:100,], Surv(days, status) ~ . )
  mod <- xgboost.ped(ped, nrounds = 100, print_every_n = 50)
  # S3 extension for stats::predict
  pred <- predict(mod, newdata = tumor[101:103,])
  expect_true(all(pred > 0))
  expect_equal(round(pred * 100, 2), c(0.69, 0.03, 2.46))
  pred_sp <- predict(mod, newdata = tumor[101:103, ], type = "surv_prob")
  expect_true(all(pred_sp >= 0 & pred_sp <= 1))
  expect_equal(round(pred_sp, 3), c(0, 0.998, 0))
  # S3 extension for pec::predictSurvProb
  sp_mat <- predictSurvProb(mod, tumor[101:103, ], times = c(1, 500))
  expect_true(all(sp_mat >= 0 & sp_mat <=1))
  expect_equal(round(sp_mat[,1], 1), rep(1, 3))
  expect_equal(round(sp_mat[,2], 1), c(.7, .9, .8))

})
