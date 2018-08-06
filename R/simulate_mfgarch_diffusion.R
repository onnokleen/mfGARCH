#' This function simulates a GARCH-MIDAS model where the short-term GARCH component is replaced by its diffusion limit, see Andersen (1998)
#' @param n.days number of days
#' @param mu mu
#' @param alpha alpha
#' @param beta beta
#' @param m m
#' @param theta theta
#' @param w1 w1
#' @param w2 w2
#' @param K K
#' @param psi psi
#' @param sigma.psi sigma.psi
#' @param low.freq low.freq
#' @param n.intraday n.intraday
#' @keywords simulate_mfgarch
#' @importFrom zoo rollapplyr
#' @importFrom stats rnorm
#' @importFrom stats setNames
#' @examples
#' \dontrun{simulate_mfgarch_diffusion(n.days = 200, mu = 0, alpha = 0.06, beta = 0.92, m = 0,
#' theta = 0.1, w1 = 1, w2 = 3, K = 12, psi = 0.98, sigma.psi = 0.1, low.freq = 10)}
#' @export
simulate_mfgarch_diffusion <- function(n.days, mu, alpha, beta, m, theta, w1 = 1, w2, K, psi, sigma.psi, low.freq = 1, n.intraday = 288) {

  if ((n.days %% low.freq) != 0) {
    stop("n.days is no multiple of low.freq")
  }

  theta_longterm <- theta

  n.days <- n.days + low.freq * K * 2

  n.sampling = 20 # 20 trades per 5 minutes
  delta <- n.intraday * n.sampling

  length_x <- n.days/low.freq
  x.innov <- rnorm(length_x, 0, sigma.psi)

  Z.p <- rnorm(n.days * delta)
  Z.sigma <- rnorm(n.days * delta)

  theta <- -log(alpha + beta)
  omega <- 1
  lambda <-
    2 * log(alpha + beta)^2 /
    (
      (((1 - (alpha + beta)^2) * (1 - beta)^2) / (alpha * (1 - beta * (alpha + beta)))) +
        6 * log(alpha + beta) +
        2 * log(alpha + beta)^2 +
        4 * (1- alpha - beta)
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
                       half.hour = rep(c(1:(n.days * n.intraday / 6)), each = 6),
                       ret = ret)

  half.hour.help <- aggregate(df.ret[c("ret")], by = list(half.hour = df.ret$half.hour, days = df.ret$days), FUN = sum)
  half.hour.vol <- aggregate(half.hour.help[c("ret")], by = list(days = half.hour.help$days), FUN = function(x) sum({x^2}))
  rm(half.hour.help)

  colnames(half.hour.vol) <- c("days", "vol")

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
