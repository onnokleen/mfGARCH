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
#' @importFrom stats setNames
#' @importFrom stats rt
#' @export
simulate_mfgarch_diffusion <- function(n.days, mu, alpha, beta, gamma, m, theta, w1 = 1, w2, K, psi, sigma.psi, low.freq = 1, student.t = NULL) {
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

  if ((n.days %% low.freq) != 0) {
    stop("n.days is no multiple of low.freq")
  }

  n.intraday <- 288

  theta_longterm <- theta

  n.days <- n.days + low.freq * K * 2

  n.sampling = 20 # 20 trades per 5 minutes
  delta <- n.intraday * n.sampling

  length_x <- n.days/low.freq
  x.innov <- rnorm(length_x, 0, sigma.psi)

  Z.p <- rnorm(n.days * delta)
  Z.sigma <- rnorm(n.days * delta)

  theta <- -log(alpha + gamma/2 + beta)
  omega <- 1
  lambda <-
    2 * log(alpha + gamma/2 + beta)^2 /
    (
      (((1 - (alpha + gamma / 2 + beta)^2) * (1 - beta)^2) / (alpha * (1 - beta * (alpha + beta + gamma / 2)))) +
        6 * log(alpha + gamma / 2 + beta) +
        2 * log(alpha + gamma / 2 + beta)^2 +
        4 * (1- alpha - beta - gamma / 2)
    )

  h <- calculate_h_andersen(ndays = n.days, delta = delta, mu = 0,
                            omega = omega, lambda = lambda, theta = theta,
                            Z = Z.sigma, pi = pi, h0 = 0.1)
  p <- calculate_p(ndays = n.days, delta = delta, mu = 0, Zp = Z.p, h = h, p0 = 0)

  p.intraday <- p[1:n.sampling == n.sampling]
  r.intraday <- p.intraday - dplyr::lag(p.intraday)



  x <- rep(0, times = length_x)
  for (ii in 2:(length_x)) {
    x[ii] <- psi * x[ii-1] + x.innov[ii]
  }
  rm(ii)

  tau <- calculate_tau(covariate = x, m = m, theta = theta_longterm, w1 = w1, w2 = w2, K = K)[-c(1:K)]

  tau <- rep(tau, each = low.freq)

  ret <- r.intraday * sqrt(rep(tau, each = n.intraday)) + mu / n.intraday

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
                    #vol_half_hour = half.hour.vol$vol,
                    real_vol = five.vol.vol$ret,
                    real_vol_half_hour = half.hour.vol$vol)

  res$real_vol_five_minutes_5_days  <- rollapplyr(res$real_vol, width = 5, FUN = mean, na.rm = TRUE, fill = NA)
  res$real_vol_five_minutes_22_days <- rollapplyr(res$real_vol, width = 22, FUN = mean, na.rm = TRUE, fill = NA)
  res$real_vol_half_hour_5_days     <- rollapplyr(res$real_vol_half_hour, width = 5, FUN = mean, na.rm = TRUE, fill = NA)
  res$real_vol_half_hour_22_days    <- rollapplyr(res$real_vol_half_hour, width = 22, FUN = mean, na.rm = TRUE, fill = NA)

  res[(low.freq * K * 2 + 1):n.days, ]
}
