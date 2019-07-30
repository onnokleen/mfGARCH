#' This function simulates a GARCH-MIDAS model. Innovations can follow a standard normal or student-t distribution.
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
#' @param low.freq number of days per low-frequency period
#' @param n.intraday number of maximum intraday returns
#' @param student.t either NULL or degrees of freedom
#' @param corr correlation between innovations (should only be used for daily tau)
#' @keywords simulate_mfgarch
#' @importFrom zoo rollapplyr
#' @importFrom stats rnorm
#' @importFrom stats rt
#' @examples
#' simulate_mfgarch(n.days = 200, mu = 0, alpha = 0.06, beta = 0.92, gamma = 0, m = 0,
#' theta = 0.1, w1 = 1, w2 = 3, K = 12, psi = 0.98, sigma.psi = 0.1, low.freq = 10)
#' @export
simulate_mfgarch <- function(n.days, mu, alpha, beta, gamma, m, theta, w1 = 1, w2, K, psi, sigma.psi, low.freq = 1, n.intraday = 288, student.t = NULL, corr = 0) {

  # n.intraday <- 288

  if (n.intraday %% 48 != 0) {
    stop("n.intraday has to be multiple of 48 (for calculating half-hour returns).")
  }

  n.days <- n.days + low.freq * K * 3

  if ((n.days %% low.freq) != 0) {
    stop("n.days is no multiple of low.freq")
  }

  if (is.null(student.t) == TRUE) {
    innov_intraday_returns <- rnorm(n.days * n.intraday)
    sim <- simulate_r(n_days = n.days, n_intraday = n.intraday,
                      alpha = alpha,
                      beta = beta,
                      gamma = gamma,
                      Z = innov_intraday_returns,
                      h0 = 1)
  } else {
    innov_intraday_returns <- rt(n = n.days * n.intraday, df = student.t) / sqrt(student.t / (student.t - 2))
    sim <- simulate_r(n_days = n.days, n_intraday = n.intraday,
                      alpha = alpha,
                      beta = beta,
                      gamma = gamma,
                      Z = innov_intraday_returns,
                      h0 = 1)
  }

  df_innov_returns <- data.frame(days = rep(c(1:n.days), each = n.intraday), innov = innov_intraday_returns /sqrt(n.intraday))
  daily_innov_returns <- aggregate(df_innov_returns[c("innov")], by = list(days = df_innov_returns$days), FUN = sum)$innov
  length_x <- n.days/low.freq
  x.innov <- corr * daily_innov_returns + sqrt(1 - corr^2) * rnorm(length_x, 0, 1)
  x.innov <- x.innov * sigma.psi
  x <- rep(0, times = length_x)
  for (ii in 2:(length_x)) {
    x[ii] <- psi * x[ii-1] + x.innov[ii]
  }
  rm(ii)
  tau <- calculate_tau(covariate = x, m = m, theta = theta, w1 = w1, w2 = w2, K = K)[-c(1:K)]
  tau <- rep(tau, each = low.freq)

  ret <- sim$ret_intraday * sqrt(rep(tau, each = n.intraday)) + mu /n.intraday

  df.ret <- data.frame(days = rep(c(1:n.days), each = n.intraday),
                       half.hour = rep(c(1:(n.days * 48)), each = n.intraday / 48),
                       ret = ret)

  half.hour.help <- aggregate(df.ret[c("ret")], by = list(half.hour = df.ret$half.hour, days = df.ret$days), FUN = sum, na.rm = TRUE)
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
                    g = sim$h_daily,
                    #vol_half_hour = half.hour.vol$vol,
                    real_vol = five.vol.vol$ret,
                    real_vol_half_hour = half.hour.vol$vol)

  res$real_vol_five_minutes_5_days  <- rollapplyr(res$real_vol, width = 5, FUN = mean, na.rm = TRUE, fill = NA)
  res$real_vol_five_minutes_22_days <- rollapplyr(res$real_vol, width = 22, FUN = mean, na.rm = TRUE, fill = NA)
  res$real_vol_half_hour_5_days     <- rollapplyr(res$real_vol_half_hour, width = 5, FUN = mean, na.rm = TRUE, fill = NA)
  res$real_vol_half_hour_22_days    <- rollapplyr(res$real_vol_half_hour, width = 22, FUN = mean, na.rm = TRUE, fill = NA)

  res[(low.freq * K * 3 + 1):n.days, ]
}

# blub <- simulate_mfgarch_rvol_dependent(n.days = 400, mu = 0, alpha = 0.06, beta = 0.92, gamma = 0, m = 0, theta = 2, w1 = 1, w2 = 3, K = 20, n.intraday = 48)
# blub

# simulate_mfgarch_rvol_dependent(n.days = 200, mu = 0, alpha = 0.06, beta = 0.92, gamma = 0, m = 0, theta = 0.1, w1 = 1, w2 = 3, K = 12, n.intraday = 1)


#' Simulate a GARCH-MIDAS similar to Wang/Ghysels with lagged RVol as covariate
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
#' @param n.intraday number of maximum intraday returns, default 288
#' @param low.freq number of days per low frequency
#' @param rvol if TRUE, the square root of the realized variance is used as a covariate
#' @examples
#' simulate_mfgarch_rv_dependent(n.days = 2200, mu = 0, alpha = 0.06, beta = 0.92, gamma = 0, m = 0,
#'   theta = 0.1, w1 = 1, w2 = 3, K = 3, low.freq = 22)
#' @export
simulate_mfgarch_rv_dependent <- function(n.days, mu, alpha, beta, gamma, m, theta, w1 = 1, w2, K, n.intraday = 288, low.freq = 1, rvol = FALSE) {
  # low.freq <- 1
  weighting <- calculate_phi(w1, w2, K)
  # n.intraday <- 288

  if (n.intraday %% 48 != 0) {
    stop("n.intraday has to be multiple of 48 (for calculating half-hour returns).")
  }

  n.days <- n.days + low.freq * K * 10

  if ((n.days %% low.freq) != 0) {
    stop("n.days is no multiple of low.freq")
  }
  innov_intraday_returns <- rnorm(n.days * n.intraday)

  sim <- simulate_r_rv_as_dependent(n_days = n.days, n_intraday = n.intraday,
                                      alpha = alpha,
                                      beta = beta,
                                      gamma = gamma,
                                      Z = innov_intraday_returns,
                                      h0 = 1, K = K, m = m, theta = theta, weights = weighting,
                                      lowfreq = low.freq,
                                      rvol = rvol)

  ret <- sim$ret_intraday + mu /n.intraday

  df.ret <- data.frame(days = rep(c(1:n.days), each = n.intraday),
                       half.hour = rep(c(1:(n.days * 48)), each = n.intraday / 48),
                       ret = ret)

  half.hour.help <- aggregate(df.ret[c("ret")], by = list(half.hour = df.ret$half.hour, days = df.ret$days), FUN = sum, na.rm = TRUE)
  half.hour.vol <- aggregate(half.hour.help[c("ret")], by = list(days = half.hour.help$days), FUN = function(x) sum({x^2}))
  rm(half.hour.help)

  colnames(half.hour.vol) <- c("days", "vol")

  five.vol.vol <- aggregate(df.ret[c("ret")], by = list(days = df.ret$days), FUN = function(x) sum(x^2))
  daily.ret    <- aggregate(df.ret[c("ret")], by = list(days = df.ret$days), FUN = sum)

  res <- data.frame(date = c(1:n.days),
                    # return_sim = sim$ret_daily,
                    return = daily.ret$ret,
                    low_freq = rep(c(1:(n.days/low.freq)), each = low.freq),
                    tau = rep(sim$tau, each = low.freq),
                    g = sim$h_daily,
                    covariate = rep(sim$rv22, each = low.freq),
                    #vol_half_hour = half.hour.vol$vol,
                    real_vol = five.vol.vol$ret,
                    real_vol_half_hour = half.hour.vol$vol)

  res$real_vol_five_minutes_5_days  <- rollapplyr(res$real_vol, width = 5, FUN = mean, na.rm = TRUE, fill = NA)
  res$real_vol_five_minutes_22_days <- rollapplyr(res$real_vol, width = 22, FUN = mean, na.rm = TRUE, fill = NA)
  res$real_vol_half_hour_5_days     <- rollapplyr(res$real_vol_half_hour, width = 5, FUN = mean, na.rm = TRUE, fill = NA)
  res$real_vol_half_hour_22_days    <- rollapplyr(res$real_vol_half_hour, width = 22, FUN = mean, na.rm = TRUE, fill = NA)

  # res$covariate <- rep(sqrt(aggregate(res[c("return")], by = list(month = res$low_freq), FUN = function(x) mean(x^2))$return), each = low.freq)

  res[(low.freq * K * 10 + 1):n.days, ]
}
