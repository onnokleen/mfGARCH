#' This function estimates a multiplicative mixed-frequency GARCH model
#' @param data data frame containing a column named date of type 'Date'.
#' @param y name of high frequency dependent variable in df.
#' @param x covariate employed in mfGARCH.
#' @param K an integer specifying lag length K in the long-term component.
#' @param low.freq a string of the low frequency variable in the df.
#' @param var.ratio.freq specify a frequency column on which the variance ratio should be calculated.
#' @param gamma if TRUE, an asymmetric GJR GARCH is used as the short-term component. If FALSE, a simple GARCH(1,1) is employed.
#' @param weighting specifies the weighting scheme employed in the long-term component. Options are "beta.one.sided" (default)
#' @keywords fit_mfgarch
#' @export fit_mfgarch
#' @importFrom magrittr %>%
#' @importFrom numDeriv jacobian
#' @importFrom stats nlminb
#' @importFrom stats optimHess
#' @importFrom dplyr select_
#' @importFrom dplyr full_join
#' @importFrom dplyr left_join
#' @importFrom dplyr tbl_df
#' @importFrom dplyr mutate_
#' @importFrom dplyr data_frame
#' @importFrom dplyr group_by_
#' @importFrom dplyr summarise_
#' @importFrom dplyr ungroup
#' @importFrom dplyr distinct
#' @importFrom stats constrOptim
#' @importFrom stats na.exclude
#' @importFrom stats pnorm
#' @importFrom stats var
#' @importFrom stats setNames
#' @importFrom stats aggregate
#' @importFrom lazyeval interp
#' @importFrom numDeriv jacobian
#' @importFrom utils tail
#' @examples fit_mfgarch(data = dplyr::filter(df_financial, date >="1973-01-01", is.na(nfci) == FALSE),
#' y = "return", x = "nfci", low.freq = "year_week", K = 52)

fit_mfgarch <- function(data, y, x = NULL, K = NULL, low.freq = "date", var.ratio.freq = NULL, gamma = TRUE, weighting = "beta.one.sided") {

  if ("date" %in% colnames(data) == FALSE) {
    stop("No date column.")
  }
  if (inherits(data$date, 'Date') == FALSE) {
    stop("Supplied date column is not of format 'Date'.")
  }
  if (is.null(x) == FALSE && K == 0) {
    warning("You specified an external covariate x but chose K = 0 - simple GARCH is estimated (K = 0).")
  }

  if (is.null(x) == TRUE) {
    warning("No external covariate x is specified - simple GARCH is estimated (K = 0).")
    x <- "date"
    K <- 0
  }
  if (is.null(K) == TRUE) {
    warning("No K is specified - simple GARCH is estimated (K = 0).")
    x <- "date"
    K <- 0
  }
  if (K < 0 || K %% 1 != 0) {
    stop("K can't be smaller than 0 and has to be an integer.")
  }
  # if ((is.null(x) == TRUE && (is.null(K) == TRUE)) || K == 0) {
  #   K <- 0
  # }
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
  if (is.null(x) == TRUE) {
    if (sum(is.na(data[[y]]) == TRUE) > 0) {
      stop(paste0("Column ", y, "contains NAs."))
    }
  } else {
    if (sum(is.na(data[[y]]) == TRUE) > 0 | sum(is.na(data[[x]]) == TRUE) > 0) {
      stop(paste0("Either column ", y, " or column ", x, "includes NAs."))
    }
  }
  if (length(unlist(unique(data[["date"]]))) != dim(data)[1]) {
    stop("There is more than one observation per high frequency (presumably date).")
  }
  if (is.null(var.ratio.freq) == FALSE) {
    if (var.ratio.freq %in% colnames(data) == FALSE) {
      stop(paste0("There is no var.ratio.freq column with name ", var.ratio.freq, "."))
    }
  }
  # Order by high frequency variable
  data <- data %>% dplyr::arrange_("date")
  # We store date in new variable because computation on integerized date seemed to be faster
  date_backup <- data[["date"]]
  data["date"] <- as.numeric(unlist(data["date"]))

  if (is.null(var.ratio.freq) == TRUE) {
    var.ratio.freq <- low.freq
    print(paste0("No frequency specified for calculating the variance ratio - default: low.freq = ", low.freq))
  }
  mutate_call_low_freq <- lazyeval::interp(~ as.integer(a), a = as.name(low.freq))
  if (x != "date") {
    df_llh <- data[,c(y, x, low.freq)] %>%
      mutate_(.dots = setNames(list(mutate_call_low_freq), low.freq))
    rm(mutate_call_low_freq)
  }

  g_zero <- var(unlist(data[y]))
  ret <- data[[y]] %>% unlist() %>% as.numeric() #unlist(df_llh[y])

  # Parameter estimation
  if (K == 0) {
    if (gamma == TRUE) {
      lf <- function(p) {
        llh_simple(y = ret,
                       mu = p["mu"],
                       alpha = p["alpha"],
                       beta = p["beta"],
                       gamma = p["gamma"],
                       m = p["m"],
                       g_zero = g_zero)
      }
      par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04, m = 0)
      ui.opt <- rbind(c(0, -1, -1, -1/2, 0), c(0, 1, 0, 0, 0), c(0, 0, 1, 0, 0))
      ci.opt <- c(-0.99999, 0, 0)
    } else {
      lf <- function(p) {
        llh_simple(y = ret,
                       mu = p["mu"],
                       alpha = p["alpha"],
                       beta = p["beta"],
                       gamma = 0,
                       m = p["m"],
                       g_zero = g_zero)

      }
      par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, m = 0)
      ui.opt <- rbind(c(0, -1, -1, 0), c(0, 1, 0, 0), c(0, 0, 1, 0))
      ci.opt <- c(-0.99999, 0, 0)
    }
    p.e.nlminb <- constrOptim(theta = par.start,
                              f = function(theta) { sum(lf(theta)) },
                              grad = NULL, ui = ui.opt, ci = ci.opt, hessian = FALSE)
    par <- p.e.nlminb$par
    returns <- data[[y]] %>% unlist() %>% as.numeric()
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
    df.fitted$residuals <- unlist((df.fitted$returns - par["mu"]) / sqrt(df.fitted$g * df.fitted$tau))
  } else { # if K > 0 we get the covariate series
    covariate <- data %>%
      select_(quote(low.freq), quote(x)) %>%
      distinct() %>%
      select_(quote(x)) %>%
      unlist()
  }

  if (K == 1) {
    if (gamma == TRUE) {
      lf <- function(p) {
        llh_mf(df = df_llh, y = ret, x = covariate, low.freq = low.freq,
               mu = p["mu"],
               omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
               alpha = p["alpha"],
               beta = p["beta"],
               gamma = p["gamma"],
               m = p["m"],
               theta = p["theta"],
               w1 = 1, w2 = 1, g_zero = g_zero, K = K)
      }
      par_start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04, m = 0, theta = 0)
      ui_opt <- rbind(c(0, -1, -1, -1/2, 0, 0), c(0, 1, 0, 0, 0, 0), c(0, 0, 1, 0, 0, 0))
      ci_opt <- c(-0.99999, 0, 0)
    } else {
      lf <- function(p) {
        llh_mf(df = df_llh,
               y = ret, x = covariate, low.freq = low.freq,
               mu = p["mu"],
               omega = 1 - p["alpha"] - p["beta"],
               alpha = p["alpha"],
               beta = p["beta"],
               gamma = 0,
               m = p["m"],
               theta = p["theta"],
               w1 = 1, w2 = 1, g_zero = g_zero, K = K)
      }
      par_start <- c(mu = 0, alpha = 0.02, beta = 0.85, m = 0, theta = 0)
      ui_opt <- rbind(c(0, -1, -1,  0, 0), c(0, 1, 0, 0, 0), c(0, 0, 1, 0, 0))
      ci_opt <- c(-0.99999, 0, 0)
    }

    p.e.nlminb <- constrOptim(theta = par_start, f = function(theta) { sum(lf(theta)) },
                              grad = NULL, ui = ui_opt, ci = ci_opt, hessian = FALSE)
    par <- p.e.nlminb$par

    tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                            theta = par["theta"], m = par["m"], w1 = 1, w2 = 1, K = K)$tau
    tau_forecast <-
      exp(sum_tau_fcts(m = par["m"],
                  i = K + 1,
                  theta = par["theta"],
                  phivar = calculate_phi(w1 = 1, w2 = 1, K = K),
                  covariate = c(data %>%
                                  select_(quote(x), quote(low.freq)) %>%
                                  distinct() %>%
                                  select_(quote(x)) %>%
                                  unlist() %>%
                                  tail(K), NA),
                  K = K))

    returns <- unlist(data[y])

    if (gamma == TRUE) {
      g <- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
             calculate_g(omega = 1 - par["alpha"] - par["beta"] - par["gamma"]/2,
                         alpha = par["alpha"], beta = par["beta"], gamma = par["gamma"],
                         as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))), g0 = g_zero))
    } else {
      g <- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
             calculate_g(omega = 1 - par["alpha"] - par["beta"],
                         alpha = par["alpha"], beta = par["beta"], gamma = 0,
                         as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))), g0 = g_zero))
    }


    df.fitted <- data %>%
      cbind(g = g, tau = tau)
    df.fitted$residuals <- unlist((df.fitted[y] - par["mu"]) / sqrt(df.fitted$g * df.fitted$tau))

  }

  if (K > 1) {
    if (gamma == TRUE) {
      lf <- function(p) {
        llh_mf(df = df_llh,
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
               K = K)
      }
      par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04, m = 0, theta = 0, w2 = 3)
      ui.opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0), c(0, 0, 0, 0, 0, 0, 1), c(0, 1, 0, 0, 0, 0, 0), c(0, 0, 1, 0, 0, 0, 0))
      ci.opt <- c(-0.99999999, 1, 0, 0)
    } else {
      lf <- function(p) {
        llh_mf(df = df_llh,
                             y = ret,
                             x = covariate,
                             low.freq = low.freq,
                             mu = p["mu"], omega = 1 - p["alpha"] - p["beta"],
                             alpha = p["alpha"],
                             beta = p["beta"],
                             gamma = 0,
                             m = p["m"],
                             theta = p["theta"],
                             w1 = 1,
                             w2 = p["w2"],
                             g_zero = g_zero,
                             K = K)
      }
      par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, m = 0, theta = 0, w2 = 3)
      ui.opt <- rbind(c(0, -1, -1, 0, 0, 0), c(0, 0, 0,  0, 0, 1), c(0, 1, 0, 0,  0, 0), c(0, 0, 1, 0, 0, 0))
      ci.opt <- c(-0.99999999, 1, 0, 0)
    }

    p.e.nlminb <- constrOptim(theta = par.start, f = function(theta) { sum(lf(theta)) },
                              grad = NULL, ui = ui.opt, ci = ci.opt, hessian = FALSE)
    par <- p.e.nlminb$par

    tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                                     w1 = 1, w2 = par["w2"],
                                     theta = par["theta"],
                                     m = par["m"], K = K)$tau

    tau_forecast <-
      exp(sum_tau_fcts(m = par["m"],
                  i = K + 1,
                  theta = par["theta"],
                  phivar = calculate_phi(w1 = 1, w2 = par["w2"], K = K),
                  covariate = c(data %>%
                                  select_(quote(x), quote(low.freq)) %>%
                                  distinct() %>%
                                  select_(quote(x)) %>%
                                  unlist() %>%
                                  tail(K), NA),
                  K = K))

    returns <- unlist(data[y])

    if (gamma == TRUE) {
      g <- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
             calculate_g(omega = 1 - par["alpha"] - par["beta"] - par["gamma"]/2,
                         alpha = par["alpha"],
                         beta = par["beta"],
                         gamma = par["gamma"],
                         as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))),
                         g0 = g_zero))
    } else {
      g <- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
             calculate_g(omega = 1 - par["alpha"] - par["beta"],
                         alpha = par["alpha"],
                         beta = par["beta"],
                         gamma = 0,
                         as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))),
                         g0 = g_zero))
    }


    df.fitted <- data %>%
      cbind(g = g, tau = tau)

    df.fitted$residuals <- unlist((df.fitted[y] - par["mu"]) / sqrt(df.fitted$g * df.fitted$tau))

  }
  df.fitted$date <- as.Date(date_backup)
  # Standard errors --------------------------------------------------------------------------------
  hessian <- solve(-optimHess(par = par, fn = function (theta) { sum(lf(theta)) }))
  rob.std.err <- sqrt(diag(hessian %*% crossprod(jacobian(func = lf, x = par)) %*% hessian))

  # Output -----------------------------------------------------------------------------------------
  output <-
    list(par = par,
         std.err = rob.std.err,
         broom.mgarch = data_frame(term = names(par),
                                   estimate = par,
                                   rob.std.err = rob.std.err,
                                   p.value = 2 * (1 - pnorm(unlist(abs(par/rob.std.err))))),
         tau = tau,
         g = g,
         df.fitted = df.fitted,
         llh = -p.e.nlminb$value,
         bic = log(sum(!is.na(tau))) * length(par) - 2 * (-p.e.nlminb$value),
         optim = p.e.nlminb)

  # Additional output if there is a long-term component (K > 0) -------------------------------------
  if (K > 0) {
    output$variance.ratio <- 100 *
                     var(log(aggregate(df.fitted$tau, by = df.fitted[var.ratio.freq],
                                       FUN = mean)[,2]),
                         na.rm = TRUE) /
                     var(log(aggregate(df.fitted$tau * df.fitted$g, by = df.fitted[var.ratio.freq],
                                       FUN = mean)[,2]),
                         na.rm = TRUE)
    output$tau_forecast = tau_forecast
  }

  # Add class mfGARCH for employing generic functions
  class(output) <- "mfGARCH"
  output
}
