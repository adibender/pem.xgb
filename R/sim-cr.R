#' Simulate competing risks time-to-event data via piece-wise exponential distribution
#'
#' Piece-wise exponential implementation of simulation algorithm described in
#' Beyersmann et al. 2009 (<doi: 10.1002/sim.3516>).
#'
#' @inheritParams pammtools::sim_pexp
#' @importFrom Formula Formula
#' @importFrom rlang is_atomic
#' @importFrom purrr reduce
#' @examples
#' library(pammtools)
#' library(dplyr)
#' library(survival)
#' # Example from Competing Risks book for comparison
#' simul.dat.cp <- function(n, h01, h02, cens.param) {
#'   times <- rexp(n, h01 + h02)
#'   ev <- rbinom(n, size = 1, prob = h01 / (h01 + h02))
#'   ev <- ifelse(ev == 0, 2, 1)
#'   cens.time <- runif(n, cens.param[1], cens.param[2])
#'   obs.time <- pmin(times, cens.time)
#'   obs.cause <- as.numeric(times == obs.time) * ev
#'   data.frame(obs.time, obs.cause)
#' }
#' set.seed(923)
#' n <- 1000
#' dat.exo7 <- simul.dat.cp(n, 0.5, 0.9, c(0.5, 3))
#' # We compute the Nelson-Aalen estimates using survfit()
#' # That’s just to get the number of events and the risk set
#' temp01 <- survfit(Surv(obs.time, obs.cause == 1) ~ 1, dat.exo7)
#' temp02 <- survfit(Surv(obs.time, obs.cause == 2) ~ 1, dat.exo7)
#' na01 <- cumsum(temp01$n.event / temp01$n.risk)
#' na02 <- cumsum(temp02$n.event / temp02$n.risk)
#'
#' plot(temp01$time, na02, type="s", ylab="Cumulative transition hazard", xlab="time")
#' lines(temp01$time, na01, type="s", col=2)
#'
#' # create data set with variables which will affect the hazard rate
#' # (not used here, but usually more complex examples of interest)
#' df <- cbind.data.frame(x1 = runif (n, -3, 3), x2 = runif (n, 0, 6)) %>%
#'  as_tibble()
#' set.seed(24032018)
#' df <- cbind.data.frame(
#'   x1 = runif (n, -3, 3),
#'   x2 = runif (n, 0, 6))
#' # two component formula specifying cause specific hazards
#' form <- ~ log(0.5)| log(0.9)
#' sim_df <- sim_pexp_cr(form, df, seq(0, 3, by =.25)) %>%
#'  mutate(
#'   cens_time = runif(n(), 0.5, 3),
#'   status = if_else(cens_time < time, 0, 1),
#'   time = pmin(time, cens_time),
#'   type = status * type)
#' temp01_2 <- survfit(Surv(time,type == 1) ~ 1, sim_df)
#' temp02_2 <- survfit(Surv(time, type == 2) ~ 1, sim_df)
#' na01_2 <- cumsum(temp01_2$n.event / temp01_2$n.risk)
#' na02_2 <- cumsum(temp02_2$n.event / temp02_2$n.risk)
#'
#' lines(temp01_2$time, na02_2, type="s", col = 3)
#' lines(temp01_2$time, na01_2, type="s", col = 4)
#' @export
sim_pexp_cr <- function(formula, data, cut) {

    Form <- Formula(formula)
    F_rhs <- attr(Form, "rhs")
    l_rhs <- length(F_rhs)
    seq_rhs <- seq_len(l_rhs)

    data <- data %>%
      mutate(
        id     = row_number(),
        time   = max(cut),
        status = 1)

  # construct eta for time-constant part
  ped  <- split_data(
      formula = Surv(time, status)~.,
      data    = select_if (data, is_atomic),
      cut     = cut,
      id      = "id") %>%
    rename("t" = "tstart")

  # calculate cause specific hazards
  for(i in seq_rhs) {
    ped[[paste0("hazard", i)]] <-  exp(eval(F_rhs[[i]], ped))
  }
  ped[["rate"]] <- reduce(ped[paste0("hazard", seq_rhs)], `+`)

  # simulate survival times
  sim_df <- ped %>%
    group_by(id) %>%
    mutate(
      time   = pammtools:::rpexp(rate = .data$rate, t = .data$t),
      status = 1L * (.data$time <= max(cut)),
      time   = pmin(.data$time, max(cut))) %>%
    filter(.data$t < .data$time & .data$time <= .data$tend)
  sim_df$type <- apply(sim_df[paste0("hazard", seq_rhs)], 1,
    function(probs)
      sample(seq_rhs, 1, prob = probs))

  sim_df %>%
    mutate(type = ifelse(.data$status == 1, .data$type, 0)) %>%
    select(-one_of(c("t", "tend", "interval", "offset", "ped_status", "rate")))

}
