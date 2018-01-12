#' This function estimates a multiplicative mixed-frequency GARCH model
#' @param n.days number of days
#' @param mu mu
#' @param alpha alpha
#' @param beta beta
#' @param gamma gamma
#' @param m m
#' @param theta theta
#' @param w1 w1
#' @param w2 w2
#' @param K K
#' @param psi psi
#' @param sigma.psi sigma.psi
#' @param low.freq low.freq
#' @param student.t student.t
#' @keywords simulate_mfgarch
#' @importFrom magrittr %>%
#' @importFrom dplyr select_
#' @importFrom dplyr mutate
#' @importFrom dplyr data_frame
#' @importFrom dplyr group_by_
#' @importFrom dplyr summarise
#' @importFrom dplyr ungroup
#' @importFrom dplyr group_by
#' @importFrom zoo rollapplyr
#' @importFrom stats rnorm
#' @importFrom stats rt
#' @export
simulate_mfgarch <- function(n.days, mu, alpha, beta, gamma, m, theta, w1 = 1, w2, K, psi, sigma.psi, low.freq = 1, student.t = NULL) {
  # Simulate a MG time series.
  #
  # Args:
  #  garch.parameters: a vector of the GARCH component parameters, i.e. mu, omega, alpha, beta, gamma.
  #  midas.parameters: .
  #  n.intraday: number of intraday returns.
  #  covariate:
  #  innovations:
  #
  # Returns:
  #  A simulated series of a MG returns, length = length(tau.innov).
  #  A series of realized volatilities based on intraday returns, length = length(covariate).
  #
  # browser()
  n.intraday <- 288

  n.days <- n.days + low.freq * K * 2

  if ((n.days %% low.freq) != 0) {
    stop("n.days is no multiple of low.freq")
  }

  length_x <- n.days/low.freq
  x.innov <- rnorm(length_x, 0, sigma.psi)

  x <- rep(0, times = length_x)
  for (ii in 2:(length_x)) {
    x[ii] <- psi * x[ii-1] + x.innov[ii]
  }
  rm(ii)

  tau <- calculate_tau(covariate = x, m = m, theta = theta, w1 = w1, w2 = w2, K = K)[-c(1:K)]

  tau <- rep(tau, each = low.freq)

  if (is.null(student.t) == TRUE) {
    sim <- simulate_r(n_days = n.days, n_intraday = n.intraday,
                      alpha = alpha,
                      beta = beta,
                      gamma = gamma,
                      Z = rnorm(n.days * n.intraday),
                      h0 = 0.1)
  } else {
    sim <- simulate_r(n_days = n.days, n_intraday = n.intraday,
                      alpha = alpha,
                      beta = beta,
                      gamma = gamma,
                      Z = rt(n = n.days * n.intraday, df = student.t),
                      h0 = 0.1)
  }

  ret <- sim$ret_intraday * sqrt(rep(tau, each = n.intraday)) + mu /n.intraday

  df.ret <- data_frame(days = rep(c(1:n.days), each = n.intraday),
                       half.hour = rep(c(1:(n.days * 48)), each = 6),
                       ret = ret)

  half.hour.vol <-
    df.ret %>%
    group_by(half.hour, days) %>%
    summarise(ret = sum(ret)) %>%
    group_by(days) %>%
    summarise(vol = sum(ret^2))

  five.vol <-
    df.ret %>%
    group_by(days) %>%
    summarise(vol = sum(ret^2),
              ret = sum(ret)) %>%
    ungroup()

  res <- data_frame(date = c(1:n.days),
                    return = five.vol$ret,
                    covariate = rep(x, each = low.freq),
                    low_freq = rep(c(1:(n.days/low.freq)), each = low.freq),
                    tau = tau,
                    g = sim$h_daily,
                    #vol_half_hour = half.hour.vol$vol,
                    real_vol = five.vol$vol,
                    real_vol_half_hour = half.hour.vol$vol) %>%
    mutate(real_vol_5_days = rollapplyr(.$real_vol, width = 5, FUN = mean, na.rm = TRUE, fill = NA)) %>%
    mutate(real_vol_22_days = rollapplyr(.$real_vol, width = 22, FUN = mean, na.rm = TRUE, fill = NA))

  res[(low.freq * K * 2 + 1):n.days, ]
}
