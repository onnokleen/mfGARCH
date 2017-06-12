#' This function estimates a multiplicative mixed-frequency GARCH model
#' @param df data frame, should contain date column
#' @param y name of high frequency dependent variable in df.
#' @param x covariate employed in mfGARCH.
#' @param K an integer specifying lag length K in the long-term component.
#' @param low.freq a string of the low frequency variable in the df.
#' @param var.ratio.freq specify a frequency column on which the variance ratio should be calculated.
#' @param gamma if TRUE, an asymmetric GJR GARCH is used as the short-term component. If FALSE, a simple GARCH(1,1) is employed.
#' @keywords fit_mfgarch
#' @export fit_mfgarch
#' @importFrom magrittr %>%
#' @importFrom numDeriv jacobian
#' @importFrom stats nlminb
#' @importFrom stats optimHess
#' @importFrom dplyr select_
#' @importFrom dplyr full_join
#' @importFrom dplyr tbl_df
#' @importFrom dplyr mutate
#' @importFrom dplyr mutate_
#' @importFrom dplyr data_frame
#' @importFrom dplyr group_by_
#' @importFrom dplyr summarise
#' @importFrom dplyr ungroup
#' @importFrom dplyr distinct
#' @importFrom stats constrOptim
#' @importFrom stats na.exclude
#' @importFrom stats var
#' @importFrom utils tail
#' @importFrom lazyeval interp
#' @importFrom numDeriv jacobian
#' @examples fit_mfgarch(data, df_mf_financial, y = "return", x = "return", low.freq = "date", K = 0)

fit_mfgarch <- function(data, y, x = NULL, K = NULL, low.freq = "date", var.ratio.freq = NULL, gamma = TRUE) {

  if (is.null(x) == TRUE || is.null(K) == TRUE) {
    print("No x or K are specified - simple GARCH is estimated (K = 0).")
    x = y
    K = 0
  }
  if (K < 0 || K %% 1 != 0) {
    stop("K can't be smaller than 0 and has to be an integer.")
  }
  if ((is.null(x) == TRUE && (is.null(K) == TRUE)) || K == 0) {
    K <- 0
  }
  if (y %in% colnames(data) == FALSE) {
    stop(paste("There is no variable in your data frame with name ", y, "."))
  }
  if (x %in% colnames(data) == FALSE && is.null(x) != FALSE) {
    stop(paste("There is no variable in your data frame with name ", x, "."))
  }
  if (low.freq %in% colnames(data) == FALSE) {
    stop(paste("There is no low freq. variable in your data frame with name ", low.freq, "."))
  }
  if ("tau" %in% colnames(data) == TRUE) {
    stop("There may not be a column named tau - it will be part of df.fitted")
  }
  if ("g" %in% colnames(data) == TRUE) {
    stop("There may not be a column named g - it will be part of df.fitted")
  }
  if (sum(is.na(select_(data, ~get(y))) == TRUE) > 0 | sum(is.na(select_(data, ~get(x))) == TRUE) > 0) {
    stop("Either y or x include NAs")
  }
  if ("date" %in% colnames(data) == FALSE) {
    stop("No date column.")
  }
  if (length(unlist(distinct(select_(data, "date")))) != dim(data)[1]) {
    stop("There is more than one observation per high frequency (presumably date).")
  }

  # If this is not the case, we may order by low frequency variable
  data <- data %>% dplyr::arrange_("date")

  data["date"] <- as.numeric(unlist(data["date"]))

  if (is.null(var.ratio.freq) == TRUE) {
    var.ratio.freq <- low.freq
    print(paste0("No frequency specified for calculating the variance ratio - default: low.freq = ", low.freq))
  }
  mutate_call_low_freq <- lazyeval::interp(~ as.integer(a), a = as.name(low.freq))
  df_llh <- data %>%
    select_(., ~get(y), ~get(x), ~get(low.freq)) %>%
    mutate_(., .dots = setNames(list(mutate_call_low_freq), low.freq))
  rm(mutate_call_low_freq)

  g_zero <- var(unlist(data[y]))
  ret <- unlist(df_llh[y])

  # Parameter estimation
  if (K == 0) {
    if (gamma == TRUE) {
      lf <- function(p) {
        sum(llh_simple(df = df_llh,
                       y = ret,
                       mu = p["mu"],
                       omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                       alpha = p["alpha"],
                       beta = p["beta"],
                       gamma = p["gamma"],
                       m = p["m"],
                       g_zero = g_zero))
      }
      par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04, m = 0)
      ui.opt <- rbind(c(0, -1, -1, -1/2, 0), c(0, 1, 0, 0, 0), c(0, 0, 1, 0, 0))
      ci.opt <- c(-0.99999, 0, 0)
    } else {
      lf <- function(p) {
        sum(llh_simple(df = df_llh,
                                 y = ret,
                                 mu = p["mu"],
                                 omega = 1 - p["alpha"] - p["beta"],
                                 alpha = p["alpha"],
                                 beta = p["beta"],
                                 gamma = 0,
                                 m = p["m"],
                                 g_zero = g_zero))
      }
      par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, m = 0)
      ui.opt <- rbind(c(0, -1, -1, 0), c(0, 1, 0, 0), c(0, 0, 1, 0))
      ci.opt <- c(-0.99999, 0, 0)
    }
    p.e.nlminb <- constrOptim(theta = par.start, f = lf, grad = NULL, ui = ui.opt, ci = ci.opt, hessian = FALSE)
    par <- p.e.nlminb$par
    returns <- data %>% select_(~get(y)) %>% unlist() %>% as.numeric()
    tau <- rep(exp(par["m"]), times = length(returns))

    if (gamma == TRUE) {
      g <- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
             calculate_g(omega = 1 - par["alpha"] - par["beta"] -  par["gamma"]/2,
                         alpha = par["alpha"],
                         beta = par["beta"],
                         gamma = par["gamma"],
                         as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))),
                         g0 = g_zero))

    } else {
      g<- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
            calculate_g(omega = 1 - par["alpha"] - par["beta"],
                        alpha = par["alpha"],
                        beta = par["beta"],
                        gamma = 0,
                        as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))),
                        g0 = g_zero))
    }

    df.fitted <-
      data_frame(returns = returns,
                 g = g,
                 tau = tau)
    df.fitted$residuals <- (df.fitted$returns - par["mu"]) / sqrt(df.fitted$g * df.fitted$tau)
  } else { # if K > 0 we get the covariate series
    covariate <- data %>%
      select_(., ~get(low.freq), ~get(x)) %>%
      distinct() %>% select_(., ~get(x)) %>%
      unlist()
  }

  if (K == 1) {
    lf <- function(p) {
      sum(likelihood_mg_mf(df = df_llh, y = y, x = x, low.freq = low.freq,
                           mu = p["mu"],
                           omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                           alpha = p["alpha"],
                           beta = p["beta"],
                           gamma = p["gamma"],
                           m = p["m"],
                           theta = p["theta"],
                           w1 = 1, w2 = 1, g_zero = g_zero, K = K))
    }

    par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04, m = 0, theta = 0)

    ui.opt <- rbind(c(0, -1, -1, -1/2, 0, 0), c(0, 1, 0, 0, 0, 0), c(0, 0, 1, 0, 0, 0))
    ci.opt <- c(-0.99999, 0, 0)

    p.e.nlminb <- constrOptim(theta = par.start, f = lf, grad = NULL, ui = ui.opt, ci = ci.opt, hessian = FALSE)
    par <- p.e.nlminb$par

    tau <- calculate_tau_mf(df = data, x = x,
                                     low.freq = low.freq,
                                     theta = par["theta"], m = par["m"], w1 = 1, w2 = 1,
                                     K = K)$tau

    returns <- unlist(data[y])

    g <- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
           calculate_g(omega = 1 - par["alpha"] - par["beta"] - par["gamma"]/2,
                       alpha = par["alpha"], beta = par["beta"], gamma = par["gamma"],
                       as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))), g0 = g_zero))

    df.fitted <- data_frame(return = returns,
                            g = g,
                            tau = tau) %>%
      mutate(., residuals = (return - par["mu"])/sqrt(g * tau))

  }

  if (K > 1) {
    lf <- function(p) {
      sum(likelihood_mg_mf(df = df_llh,
                           y = ret,
                           x = covariate,
                           low.freq = low.freq,
                           mu = p["mu"], omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                           alpha = p["alpha"],
                           beta = p["beta"],
                           gamma = p["gamma"],
                           m = p["m"],
                           theta = p["theta"],
                           w1 = 1,
                           w2 = p["w2"],
                           g_zero = g_zero,
                           K = K))
    }

    par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04, m = 0, theta = 0, w2 = 3)

    ui.opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0), c(0, 0, 0, 0, 0, 0, 1), c(0, 1, 0, 0, 0, 0, 0), c(0, 0, 1, 0, 0, 0, 0))
    ci.opt <- c(-0.99999999, 1, 0, 0)

    p.e.nlminb <- constrOptim(theta = par.start, f = lf, grad = NULL, ui = ui.opt, ci = ci.opt, hessian = FALSE)
    par <- p.e.nlminb$par

    tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                                     w1 = 1, w2 = par["w2"],
                                     theta = par["theta"],
                                     m = par["m"], K = K)$tau

    tau_forecast <-
      exp(sum_tau(m = par["m"],
                  i = K + 1,
                  theta = par["theta"],
                  phivar = calculate_phi(w1 = 1, w2 = par["w2"], K = K),
                  covariate = c(data %>% select_(., ~get(x), ~get(low.freq)) %>%
                                  distinct() %>% select_(~get(x)) %>%
                                  unlist() %>%
                                  tail(K), NA),
                  K = K))

    returns <- unlist(data[y])

    g <- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
                       calculate_g(omega = 1 - par["alpha"] - par["beta"] - par["gamma"]/2,
                                   alpha = par["alpha"],
                                   beta = par["beta"],
                                   gamma = par["gamma"],
                                   as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))),
                                   g0 = g_zero))

    df.fitted <- data %>%
      cbind(g = g, tau = tau) %>%
      mutate(., residuals = (returns - par["mu"])/sqrt(g * tau))

  }

  if (K == 0) {
    if (gamma == TRUE) {
      lf_mgarch_score <- function(p) {
        llh_simple(df = df_llh,
                             y = ret,
                             mu = p["mu"],
                             omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                             alpha = p["alpha"],
                             beta = p["beta"],
                             gamma = p["gamma"],
                             m = p["m"],
                             g_zero = g_zero)
      }
    } else {
      lf_mgarch_score <- function(p) {
        llh_simple(df = df_llh,
                             y = ret,
                             mu = p["mu"],
                             omega = 1 - p["alpha"] - p["beta"],
                             alpha = p["alpha"],
                             beta = p["beta"],
                             gamma = 0,
                             m = p["m"],
                             g_zero = g_zero)
      }
    }
  }

  if (K == 1) {
    lf_mgarch_score <- function(p) {
      likelihood_mg_mf(df = df_llh,
                       y = y,
                       x = x,
                       low.freq = low.freq,
                       mu = p["mu"], omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                       alpha = p["alpha"],
                       beta = p["beta"],
                       gamma = p["gamma"],
                       m = p["m"],
                       theta = p["theta"],
                       w1 = 1, w2 = 1, g_zero = g_zero, K = K)
    }
  }

  if (K > 1) {
    lf_mgarch_score <- function(p) {
      likelihood_mg_mf(df = df_llh,
                       y = ret,
                       x = covariate,
                       low.freq = low.freq,
                       mu = p["mu"],
                       omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                       alpha = p["alpha"],
                       beta = p["beta"],
                       gamma = p["gamma"],
                       m = p["m"],
                       theta = p["theta"],
                       w1 = 1, w2 = p["w2"],
                       g_zero = g_zero, K = K)
    }
  }

  hessian <- solve(-optimHess(par = par, fn = lf))
  G <- jacobian(func = lf_mgarch_score, x = par)
  GGsum <- t(G) %*% G
  rob.std.err <- sqrt(diag(hessian %*% GGsum %*% hessian))

  if (K == 0) {
    if (gamma == TRUE) {
      output <-
        list(par = par,
             std.err = rob.std.err[1:5],
             broom.mgarch = data_frame(term = c("mu", "alpha", "beta", "gamma", "m"),
                                       estimate = c(par[1:5]),
                                       rob.std.err = c(rob.std.err[1:5])),
             g = g,
             tau = tau,
             df.fitted = df.fitted,
             llh = -p.e.nlminb$value,
             bic = log(sum(!is.na(tau))) * 5 - 2 * (-p.e.nlminb$value),
             optim = p.e.nlminb)
    } else {
      output <-
        list(par = par,
             std.err = rob.std.err[1:4],
             broom.mgarch = data_frame(term = c("mu", "alpha", "beta", "m"),
                                       estimate = c(par[1:4]),
                                       rob.std.err = c(rob.std.err[1:4])),
             g = g,
             tau = tau,
             df.fitted = df.fitted,
             llh = -p.e.nlminb$value,
             bic = log(sum(!is.na(tau))) * 4 - 2 * (-p.e.nlminb$value),
             optim = p.e.nlminb)
    }
  }

  if (K == 1) {
    output <-
      list(mgarch = p.e.nlminb,
           par = par,
           std.err = rob.std.err[1:6],
           broom.mgarch = data_frame(term = c("mu", "alpha", "beta", "gamma", "m", "theta"),
                                     estimate = par[1:6],
                                     rob.std.err = rob.std.err[1:6]),
           tau = tau,
           g = g,
           df.fitted = df.fitted,
           llh = -p.e.nlminb$value,
           bic = log(sum(!is.na(tau))) * 6 - 2 * (-p.e.nlminb$value))
  }

  if (K > 1) {
    output <-
      list(mgarch = p.e.nlminb,
           broom.mgarch = data_frame(term = c("mu", "alpha", "beta", "gamma", "m", "theta", "w2"),
                                     estimate = c(par[1:7]),
                                     rob.std.err = c(rob.std.err[1:7])),
           par = par,
           std.err = rob.std.err[1:7],
           tau = tau,
           g = g,
           df.fitted = df.fitted,
           variance.ratio = 100 * (df.fitted %>%
                                     group_by_(var.ratio.freq) %>%
                                     summarise(mean.tau = mean(tau), mean.tau.g = mean(tau * g)) %>%
                                     ungroup() %>%
                                     summarise(var.ratio = var(log(mean.tau), na.rm = TRUE)/var(log(mean.tau.g), na.rm = TRUE)) %>%
                                     unlist()),
           tau_forecast = tau_forecast,
           llh = -p.e.nlminb$value,
           bic = log(sum(!is.na(tau))) * 7 - 2 * (-p.e.nlminb$value))
  }

  class(output) <- "mfGARCH"
  output
}


#' @keywords internal
calculate_tau_mf <- function(df, x, low.freq, w1, w2, theta, m, K) {
  phi.var <- calculate_phi(w1, w2, K)
  covariate <- c(rep(NA, times = K), x)
  tau <- c(rep(NA, times = K), exp(sum_tau_zwei(m = m, theta = theta, phivar = phi.var, covariate = x, K = K)))
  left_join(df, cbind(unique(df[low.freq]), tau), by = low.freq)
}

#' @keywords internal
likelihood_mg_mf <- function(df, x, y, low.freq, mu, omega, alpha, beta, gamma, m, theta, w1 = 1, w2 = 1, g_zero, K = 2) {

  tau <- calculate_tau_mf(df = df, x = x, low.freq = low.freq, w1 = w1, w2 = w2, theta = theta, m = m, K = K)$tau
  ret <- y
  ret <- ret[which.min(is.na(tau)):length(ret)]  # lags can't be used for likelihood
  tau <- tau[which.min(is.na(tau)):length(tau)]
  g <- calculate_g(omega = omega, alpha = alpha, beta = beta, gamma = gamma,
                   returns = ((ret - mu)/sqrt(tau)), g0 = g_zero)

  if (sum(g <= 0) > 0) {
    rep(NA, times = length(y))
  } else {
    1/2 * log(2 * pi) + 1/2 * log(g * tau) + 1/2 * (ret - mu)^2/(g * tau)
  }
}

#' @keywords internal
llh_simple <- function(df, y, mu, omega, alpha, beta, gamma, m, g_zero) {
  ret <- y
  ret_std <- (ret - mu)/sqrt(exp(m))
  g <- calculate_g(omega = omega, alpha = alpha, beta = beta, gamma = gamma,
                   returns = ret_std, g0 = g_zero)
  1/2 * log(2 * pi) + 1/2 * log(g * exp(m)) + 1/2 * (ret - mu)^2/(g * exp(m))
}
