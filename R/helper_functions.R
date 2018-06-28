#' @keywords internal
forecast_garch <- function(omega, alpha, beta, gamma, g, ret, steps.ahead) {
  omega / (1 - alpha - gamma/2 - beta) + (alpha + beta + gamma/2)^(steps.ahead - 1) * (omega + (alpha + gamma/2 * as.numeric(ret < 0)) * ret^2 + beta * g - omega / (1 - alpha - gamma/2 - beta))
}


# Calculates the phi weighting scheme
#' @keywords internal
calculate_phi <- function(w1, w2, K) {
    weights <- sapply(c(1:K),
                      FUN = function(j) (j / (K + 1))^(w1 - 1) * (1 - j / (K + 1))^(w2 - 1))
    weights <- weights/sum(weights)
    weights
}


#' @keywords internal
calculate_tau <- function(covariate, w1, w2, theta, m, K) { # used for simulation
    phi_var <- calculate_phi(w1, w2, K)
    covariate <- c(rep(NA, times = K), covariate)
    tau <- c(rep(NA, times = K),
             exp(sum_tau(m = m, theta = theta, phivar = phi_var, covariate = covariate, K = K)))
    tau
}

#' @keywords internal
calculate_tau_mf <- function(df, x, low.freq, w1, w2, theta, m, K,
                             x.two = NULL, K.two = NULL, theta.two = NULL,
                             low.freq.two = NULL, w1.two = NULL, w2.two = NULL) {
  phi.var <- calculate_phi(w1, w2, K)
  covariate <- c(rep(NA, times = K), x)
  tau <- c(rep(NA, times = K),
           exp(sum_tau(m = m, theta = theta, phivar = phi.var, covariate = x, K = K)))

  result <- merge(df, cbind(unique(df[low.freq]), tau), by = low.freq)

  if (is.null(x.two) == FALSE) {
    phi.var.two <- calculate_phi(w1.two, w2.two, K.two)
    covariate.two <- c(rep(NA, times = K.two), x.two)
    tau.two <- c(rep(NA, times = K.two),
                 exp(sum_tau(m = 0, theta = theta.two, phivar = phi.var.two,
                             covariate = x.two, K = K.two)))
    result <- merge(result, cbind(unique(df[low.freq.two]), tau.two), by = low.freq.two)

    result$tau.one <- result$tau # store tau component due to first covariate
    result$tau <- result$tau.one * result$tau.two # generate joint tau component
  }

  result
}

# df_test <- df_mfgarch %>% filter(is.na(vix) == FALSE)
# llh_mf(df_test, y = df_test$return, x = unlist(unique(df_test[c("date", "vix")])["vix"]), K = 3, mu = 0, m = 0, theta = 0.1, low.freq = "date", omega = 0.01, alpha = 0.06, beta = 0.9, gamma = 0, g_zero = 1, x.two = unlist(unique(df_test[c("year_month", "dindpro")])["dindpro"]), K.two = 12, theta.two = 0.1, low.freq.two = "year_month", w1.two = 1, w2.two = 1)
#' @keywords internal
llh_mf <-
  function(df, x, y, low.freq, mu, omega, alpha, beta, gamma,
           m, theta, w1 = 1, w2 = 1, g_zero, K,
           x.two = NULL, K.two = NULL, theta.two = NULL,
           low.freq.two = NULL, w1.two = NULL, w2.two = NULL) {

    if (is.null(x.two) == FALSE) {
      tau <- calculate_tau_mf(df = df, x = x, low.freq = low.freq,
                              w1 = w1, w2 = w2, theta = theta, m = m, K = K,
                              x.two = x.two, K.two = K.two, theta.two = theta.two,
                              low.freq.two = low.freq.two,
                              w1.two = w1.two, w2.two = w2.two)$tau
    } else {
      tau <- calculate_tau_mf(df = df, x = x, low.freq = low.freq,
                              w1 = w1, w2 = w2, theta = theta, m = m, K = K)$tau
    }

    ret <- y
    ret <- ret[which.min(is.na(tau)):length(ret)]  # lags can't be used for likelihood
    tau <- tau[which.min(is.na(tau)):length(tau)]
    g <- calculate_g(omega = omega, alpha = alpha, beta = beta, gamma = gamma,
                     returns = ((ret - mu)/sqrt(tau)), g0 = g_zero)

    if (sum(g <= 0) > 0) {
      #rep(NA, times = length(y))
      #stop("g_t seems to be negative for at least one point in time?")
      rep(NA, times = length(g))
    } else {
      1/2 * log(2 * pi) + 1/2 * log(g * tau) + 1/2 * (ret - mu)^2/(g * tau)
    }
  }

#' @keywords internal
llh_simple <- function(y, mu, alpha, beta, gamma, m, g_zero) {
  omega <- 1 - alpha - beta - gamma / 2
  ret <- y
  ret_std <- (ret - mu)/sqrt(exp(m))
  g <- calculate_g(omega = omega, alpha = alpha, beta = beta, gamma = gamma,
                   returns = ret_std, g0 = g_zero)
  1/2 * log(2 * pi) + 1/2 * log(g * exp(m)) + 1/2 * (ret - mu)^2/(g * exp(m))
}

