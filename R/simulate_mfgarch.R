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
#' @importFrom zoo rollapplyr
#' @importFrom stats rnorm
#' @importFrom stats rt
#' @example simulate_mfgarch(n.days = 2000, mu = 0, alpha = 0.06, beta = 0.92, gamma = 0, m = 0,
#' theta = 0.1, w1 = 1, w2 = 3, K = 12, psi = 0.98, sigma.psi = 0.1, low.freq = 100, student.t = NULL)
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

  df.ret <- data.frame(days = rep(c(1:n.days), each = n.intraday),
                       half.hour = rep(c(1:(n.days * 48)), each = 6),
                       ret = ret)

  half.hour.help <- aggregate(df.ret[c("ret")], by = list(half.hour = df.ret$half.hour, days = df.ret$days), FUN = sum)
  half.hour.vol <- aggregate(half.hour.help[c("ret")], by = list(days = half.hour.help$days), FUN = function(x) sum({x^2}))
  rm(half.hour.help)

  colnames(half.hour.vol) <- c("days", "vol")
  #sum(half.hour.vol != half.hour.vol.2, na.rm = TRUE)

  # deprecated
  # half.hour.vol <-
  #   df.ret %>%
  #   group_by(half.hour, days) %>%
  #   summarise(ret = sum(ret)) %>%
  #   group_by(days) %>%
  #   summarise(vol = sum(ret^2))

  five.vol.vol <- aggregate(df.ret[c("ret")], by = list(days = df.ret$days), FUN = function(x) sum(x^2))
  daily.ret    <- aggregate(df.ret[c("ret")], by = list(days = df.ret$days), FUN = sum)

  res <- data.frame(date = c(1:n.days),
                    return = daily.ret$ret,
                    covariate = rep(x, each = low.freq),
                    low_freq = rep(c(1:(n.days/low.freq)), each = low.freq),
                    tau = tau,
                    g = sim$h_daily,
                    #vol_half_hour = half.hour.vol$vol,
                    real_vol = five.vol.vol$ret,
                    real_vol_half_hour = half.hour.vol$vol)

  res$real_vol_five_minutes_5_days  <- rollapplyr(res$real_vol, width = 5, FUN = mean, na.rm = TRUE, fill = NA)
  res$real_vol_five_minutes_22_days <- rollapplyr(res$real_vol, width = 22, FUN = mean, na.rm = TRUE, fill = NA)
  res$real_vol_half_hour_5_days     <- rollapplyr(res$real_vol_half_hour, width = 5, FUN = mean, na.rm = TRUE, fill = NA)
  res$real_vol_half_hour_22_days    <- rollapplyr(res$real_vol_half_hour, width = 22, FUN = mean, na.rm = TRUE, fill = NA)

  res[(low.freq * K * 2 + 1):n.days, ]
}
