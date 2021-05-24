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


test_that("competing risks setting words", {

  data("cancer", package = "survival")
  mgus2$time <- with(mgus2, ifelse(pstat == 0, futime, ptime))
  mgus2$status <- with(mgus2, ifelse(pstat == 0, 2 * death, 1))
  mgus2 <- mgus2 %>%
    dplyr::select(-id, -ptime, -futime, -death, -pstat) %>%
    dplyr::mutate_if(is.numeric, ~ifelse(is.na(.), mean(., na.rm = TRUE), .))

  ped <- as_ped(mgus2, Surv(time, status)~.)
  xgb1 <- xgb.train.ped(
    params  = list(eta = .3),
    data    = dplyr::filter(ped, cause == 1),
    nrounds = 20L,
    verbose = FALSE)
  xgb2 <- xgb.train.ped(
    params  = list(eta = .3),
    data    = dplyr::filter(ped, cause == 2),
    nrounds = 20L,
    verbose = FALSE)
  xgb_list <- list(xgb1, xgb2)
  event_prob1 <- predictEventProb(
    object = xgb_list,
    newdata = mgus2[1:10,],
    times = c(0.1, 3, 5, 200),
    cause = 1)
  event_prob2 <- predictEventProb(
    object = xgb_list,
    newdata = mgus2[1:10,],
    times = c(0.1, 3, 5, 200),
    cause = 2)


})
