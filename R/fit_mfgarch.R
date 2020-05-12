#' This function estimates a multiplicative mixed-frequency GARCH model. For the sake of numerical stability, it is best to multiply log returns by 100.
#' @param data data frame containing a column named date of type 'Date'.
#' @param y name of high frequency dependent variable in df.
#' @param x covariate employed in mfGARCH.
#' @param K an integer specifying lag length K in the long-term component.
#' @param low.freq a string of the low frequency variable in the df.
#' @param var.ratio.freq specify a frequency column on which the variance ratio should be calculated.
#' @param gamma if TRUE, an asymmetric GJR-GARCH is used as the short-term component. If FALSE, a simple GARCH(1,1) is employed.
#' @param weighting specifies the weighting scheme employed in the long-term component. Options are "beta.restricted" (default) or "beta.unrestricted"
#' @param x.two optional second covariate
#' @param K.two lag lgenth of optional second covariate
#' @param low.freq.two low frequency of optional second covariate
#' @param weighting.two specifies the weighting scheme employed in the optional second long-term component. Currently, the only option is "beta.restricted"
#' @param multi.start if TRUE, optimization is carried out with multiple starting values
#' @param control a list
#' @keywords fit_mfgarch
#' @export
#'
#' @return A list of class mfGARCH with letters and numbers.
#' \itemize{
#' \item par - vector of estimated parameters
#' \item rob.std.err - sandwich/HAC-type standard errors
#' \item broom.mgarch - a broom-like data.frame with entries
#' 1) estimate: column of estimated parameters
#' 2) rob.std.err - sandwich/HAC-type standard errors
#' 3) p.value - p-values derived from sandwich/HAC-type standard errors
#' 4) opg.std.err - Bollerslev-Wooldrige/OPG standard errors for GARCH processes
#' 5) opg.p.value - corresponding alternative p-values
#' \item tau - fitted long-term component
#' \item g - fitted short-term component
#' \item df.fitted - data frame with fitted values and residuals
#' \item K - chosen lag-length in the long-term component
#' \item weighting.scheme - chosen weighting scheme
#' \item llh - log-likelihood value at estimated parameter vector
#' \item bic - corresponding BIC value
#' \item y - dependent variable y
#' \item optim - output of the optimization routine
#' \item K.two - lag-lenth of x.two if two covariates are employed
#' \item weighting.scheme.two - chosen weighting scheme of x.two (if K.two != NULL)
#' \item tau.forecast - one-step ahead forecast of the long-term component
#' \item variance.ratio - calculated variance ratio
#' \item est.weighting - estimated weighting scheme
#' \item est.weighting.two - estimated weighting scheme of x.two (if K.two != NULL)
#' }
#'
#' @importFrom numDeriv jacobian
#' @importFrom stats nlminb
#' @importFrom numDeriv hessian
#' @importFrom stats constrOptim
#' @importFrom stats na.exclude
#' @importFrom stats optim
#' @importFrom stats pnorm
#' @importFrom stats var
#' @importFrom stats aggregate
#' @importFrom numDeriv jacobian
#' @importFrom maxLik maxLik
#' @importFrom utils tail
#' @examples
#' \dontrun{
#' fit_mfgarch(data = df_financial, y = "return", x = "nfci", low.freq = "week", K = 52)
#' fit_mfgarch(data = df_mfgarch, y = "return", x = "nfci", low.freq = "year_week", K = 52,
#' x.two = "dindpro", K.two = 12, low.freq.two = "year_month", weighting.two = "beta.restricted")
#' }

fit_mfgarch <- function(data, y, x = NULL, K = NULL, low.freq = "date", var.ratio.freq = NULL, gamma = TRUE, weighting = "beta.restricted", x.two = NULL, K.two = NULL, low.freq.two = NULL, weighting.two = NULL, multi.start = FALSE, control = list(par.start = NULL)) {

  print("For ensuring numerical stability of the parameter optimization and inversion of the Hessian, it is best to multiply log returns by 100.")

  if (is.null(weighting.two) == FALSE) {
    if (weighting.two != "beta.restricted") {
      stop("Right now, only beta.restricted weighting scheme is employed for the second covariate.")
    }
  }

  if (is.null(x.two) == FALSE) {
    weighting.two <- "beta.restricted"
  }

  if (is.null(x.two) == FALSE && gamma == FALSE) {
    stop("Regarding two covariates, only asymmetric GJR-GARCH component is implemented.")
  }

  # if (K == 1 && is.null(x.two) == FALSE && K.two != 1) {
  #   stop("Regarding two covariates, only K.two = 1 is implemented if K = 1.")
  # }
  if (is.null(K.two) == FALSE) {
    if (K == 1 & K.two > 1) {
      stop("Regarding two covariates with one of them being equal to one, only K.two = 1 is implemented.")
    }
  }

  if (is.null(x.two) == FALSE) {
    print("Specifying two covariates may lead to long estimation times.")
  }

  if (weighting %in% c("beta.restricted", "beta.unrestricted") == FALSE) {
    stop("Incorrect weighting scheme specified - options are \"beta.restricted\" and \"beta.unrestricted\".")
  }
  if (gamma %in% c(TRUE, FALSE) == FALSE) {
    stop("Gamma can't be anything different than TRUE or FALSE.")
  }
  if ("date" %in% colnames(data) == FALSE) {
    stop("No date column.")
  }
  if (inherits(data$date, 'Date') == FALSE) {
    stop(paste0("Supplied date column is not of format 'Date'. It is of class '", class(data$date), "'."))
  }
  if (inherits(data[[low.freq]], 'Date') == FALSE) {
    stop(paste0("Supplied low.freq column is not of format 'Date'. It is of class '", class(data[[low.freq]]), "'."))
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
  if (dim(unique(data[c(x, low.freq)]))[1] > dim(unique(data[c(low.freq)]))[1]) {
    stop("There is more than one unique observation per low frequency entry.")
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
  data <- data[order(data$date), ]
  # Deprecated dplyr version
  #data <- dplyr::arrange_(data, "date")
  # We store date in new variable because computation on integerized date seemed to be faster
  date_backup <- data[["date"]]
  data["date"] <- as.numeric(unlist(data["date"]))

  if (is.null(var.ratio.freq) == TRUE) {
    var.ratio.freq <- low.freq
    print(paste0("No frequency specified for calculating the variance ratio - default: low.freq = ", low.freq))
  }

  low_freq_backup <- data[, low.freq]
  if (x != "date") {
    if (is.null(x.two) == TRUE) {
      df_llh <- data[, c(y, x, low.freq)]
      df_llh[, low.freq] <- as.integer(unlist(df_llh[ , low.freq]))
    } else {
      low_freq.two_backup <- data[, low.freq.two]
      if (low.freq != low.freq.two) { # if they are different, both have to be included in df_llh
        df_llh <- data[, c(y, x, low.freq, x.two, low.freq.two)]
        df_llh[, low.freq] <- as.integer(unlist(df_llh[ , low.freq]))
        df_llh[, low.freq.two] <- as.integer(unlist(df_llh[ , low.freq.two]))
      } else { # else, the low.freq column is needed only once
        df_llh <- data[, c(y, x, low.freq, x.two)]
        df_llh[, low.freq] <- as.integer(unlist(df_llh[ , low.freq]))
      }
    }
  }

  g_zero <- var(unlist(data[[y]]))
  ret <- data[[y]]

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

    if(is.null(control$par.start) == FALSE) {
      par.start <- control$par.start
    }

    p.e.nlminb <- constrOptim(theta = par.start,
                              f = function(theta) { sum(lf(theta)) },
                              grad = NULL, ui = ui.opt, ci = ci.opt, hessian = FALSE)

    if (multi.start == TRUE && gamma == TRUE) {
      p.e.nlminb.two <- try({
        suppressWarnings(optim(par = p.e.nlminb$par, fn = function (theta) {
          if( is.na(sum(lf(theta))) == TRUE) {
            NA
          } else {
            sum(lf(theta))
          }
        }, method = "BFGS"))}, silent = TRUE)

      if (class(p.e.nlminb.two) == "try-error") {
        print("Second-step BFGS optimization failed. Fallback: First-stage Nelder-Mead estimate.")
      } else {
        if (p.e.nlminb.two$value < p.e.nlminb$value) {
          p.e.nlminb <- p.e.nlminb.two
        }
      }
    }

    p.e.nlminb$value <- -p.e.nlminb$value

    par <- p.e.nlminb$par
    returns <- as.numeric(unlist(data[[y]]))
    tau <- rep(exp(par["m"]), times = length(returns))

    if (gamma == TRUE) {
      g <- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
             calculate_g(omega = 1 - par["alpha"] - par["beta"] -  par["gamma"]/2,
                         alpha = par["alpha"],
                         beta = par["beta"],
                         gamma = par["gamma"],
                         as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))),
                         g0 = g_zero))
      tau <- rep(exp(par["m"]), times = length(g))

    } else {
      g <- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
            calculate_g(omega = 1 - par["alpha"] - par["beta"],
                        alpha = par["alpha"],
                        beta = par["beta"],
                        gamma = 0,
                        as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))),
                        g0 = g_zero))
      tau <- rep(exp(par["m"]), times = length(g))
    }

    if ((var.ratio.freq %in% c("date", "low.freq")) == FALSE) {
      df.fitted <- cbind(data[c("date", y, var.ratio.freq)], g = g, tau = tau)
    } else {
      df.fitted <- cbind(data[c("date", y)], g = g, tau = tau)
    }

    df.fitted$residuals <- unlist((df.fitted[y] - par["mu"]) / sqrt(df.fitted$g * df.fitted$tau))
  } else { # if K > 0 we get the covariate series
    covariate <- unlist(unique(data[c(low.freq, x)])[x])

    if (is.null(x.two) == FALSE) {
      covariate.two <- unlist(unique(data[c(low.freq.two, x.two)])[x.two])
    }
  }

  if (K == 1) {
    if (is.null(K.two) == FALSE) {
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
                 w1 = 1, w2 = 1, g_zero = g_zero, K = K,
                 x.two = covariate.two,
                 K.two = K.two, low.freq.two = low.freq.two,
                 theta.two = p["theta.two"], w1.two = 1, w2.two = 1)
        }
        par_start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04, m = 0, theta = 0, theta.two = 0)
        ui_opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0), c(0, 1, 0, 0, 0, 0, 0), c(0, 0, 1, 0, 0, 0, 0))
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
                 w1 = 1, w2 = 1, g_zero = g_zero, K = K,
                 x.two = covariate.two,
                 K.two = K.two, low.freq.two = low.freq.two,
                 theta.two = p["theta.two"], w1.two = 1, w2.two = 1)
        }
        par_start <- c(mu = 0, alpha = 0.02, beta = 0.85, m = 0, theta = 0, theta.two = 0)
        ui_opt <- rbind(c(0, -1, -1,  0, 0, 0), c(0, 1, 0, 0, 0, 0), c(0, 0, 1, 0, 0, 0))
        ci_opt <- c(-0.99999, 0, 0)
      }

    } else {

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
    }

    if(is.null(control$par.start) == FALSE) {
      par.start <- control$par.start
    }

    p.e.nlminb <- constrOptim(theta = par_start, f = function(theta) { sum(lf(theta)) },
                              grad = NULL, ui = ui_opt, ci = ci_opt, hessian = FALSE)
    par <- p.e.nlminb$par
    p.e.nlminb$value <- -p.e.nlminb$value

    if (is.null(x.two) == FALSE) {
      tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                              w1 = 1, w2 = 1, theta = par["theta"], m = par["m"], K = K,
                              x.two = covariate.two, K.two = K.two, theta.two = par["theta.two"],
                              low.freq.two = low.freq.two,
                              w1.two = 1, w2.two = 1)$tau
    } else {
      tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                              theta = par["theta"], m = par["m"], w1 = 1, w2 = 1, K = K)$tau
    }

    tau_forecast <-
      exp(sum_tau_fcts(m = par["m"],
                       i = K + 1,
                       theta = par["theta"],
                       phivar = calculate_phi(w1 = 1, w2 = 1, K = K),
                       covariate = c(tail(unlist(unique(data[c(x, low.freq)])[x]), K), NA),
                       K = K))

    if (is.null(x.two) == FALSE) {
      tau_forecast <-
        tau_forecast *
        exp(sum_tau_fcts(m = 0,
                         i = K.two + 1,
                         theta = par["theta.two"],
                         phivar = calculate_phi(w1 = 1, w2 = 1, K = K.two),
                         covariate = c(tail(unlist(unique(data[c(x.two, low.freq.two)])[x.two]), K.two), NA),
                         K = K.two))
    }


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

    if ((var.ratio.freq %in% c("date", "low.freq")) == FALSE) {
      df.fitted <- cbind(data[c("date", y, low.freq, x, var.ratio.freq)], g = g, tau = tau)
    } else {
      df.fitted <- cbind(data[c("date", y, low.freq, x)], g = g, tau = tau)
    }

    df.fitted$residuals <- unlist((df.fitted[y] - par["mu"]) / sqrt(df.fitted$g * df.fitted$tau))

  }

  if (K > 1) {
    if (gamma == TRUE) {
      if (weighting == "beta.restricted" & is.null(K.two) == TRUE) {
        lf <- function(p) {
          llh_mf(df = df_llh,
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
                 w1 = 1,
                 w2 = p["w2"],
                 g_zero = g_zero,
                 K = K)
        }
        par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04,
                       m = 0, theta = 0, w2 = 3)
        ui.opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0),
                        c(0,  0,  0,    0, 0, 0, 1),
                        c(0,  1,  0,    0, 0, 0, 0),
                        c(0,  0,  1,    0, 0, 0, 0))
        ci.opt <- c(-0.99999999, 1, 0, 0)
      }
      if (weighting == "beta.restricted" & is.null(K.two) == FALSE) {
        if (K.two == 1) {
          lf <- function(p) {
            llh_mf(df = df_llh,
                   y = ret,
                   x = covariate,
                   low.freq = low.freq,
                   mu = p["mu"],
                   omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                   alpha = p["alpha"], beta = p["beta"], gamma = p["gamma"],
                   m = p["m"], theta = p["theta"],
                   w1 = 1, w2 = p["w2"], g_zero = g_zero, K = K,
                   x.two = covariate.two,
                   K.two = 1, low.freq.two = low.freq.two,
                   theta.two = p["theta.two"], w1.two = 1, w2.two = 1)
          }
          par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04,
                         m = 0, theta = 0, w2 = 3, theta.two = 0)
          ui.opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0, 0),
                          c(0,  0,  0,    0, 0, 0, 1, 0),
                          c(0,  1,  0,    0, 0, 0, 0, 0),
                          c(0,  0,  1,    0, 0, 0, 0, 0))
          ci.opt <- c(-0.99999999, 1, 0, 0)
        }
        if (K.two > 1) {
          if (weighting.two == "beta.restricted") {
            lf <- function(p) {
              llh_mf(df = df_llh,
                     y = ret,
                     x = covariate,
                     low.freq = low.freq,
                     mu = p["mu"],
                     omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                     alpha = p["alpha"], beta = p["beta"], gamma = p["gamma"],
                     m = p["m"], theta = p["theta"],
                     w1 = 1, w2 = p["w2"], g_zero = g_zero, K = K,
                     x.two = covariate.two,
                     K.two = K.two, low.freq.two = low.freq.two,
                     theta.two = p["theta.two"], w1.two = 1, w2.two = p["w2.two"])
            }
            par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04,
                           m = 0, theta = 0, w2 = 3, theta.two = 0, w2.two = 3)
            ui.opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0, 0, 0),
                            c(0,  0,  0,    0, 0, 0, 1, 0, 0),
                            c(0,  1,  0,    0, 0, 0, 0, 0, 0),
                            c(0,  0,  1,    0, 0, 0, 0, 0, 0),
                            c(0,  0,  0,    0, 0, 0, 0, 0, 1))
            ci.opt <- c(-0.99999999, 1, 0, 0, 1)
          }
          if (weighting.two != "beta.restricted") {
            stop("Weighting scheme for second variable can only be beta.restricted.")
          }
        }
      }

      if (weighting == "beta.unrestricted" & is.null(K.two) == TRUE){
        lf <- function(p) {
          llh_mf(df = df_llh,
                 y = ret,
                 x = covariate,
                 low.freq = low.freq,
                 mu = p["mu"],
                 omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                 alpha = p["alpha"], beta = p["beta"], gamma = p["gamma"],
                 m = p["m"], theta = p["theta"], w1 = p["w1"], w2 = p["w2"],
                 g_zero = g_zero,
                 K = K)
        }
        par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04,
                       m = 0, theta = 0, w1 = 1.0000001, w2 = 3)
        ui.opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0, 0),
                        c(0,  0,  0,  0,   0, 0, 1, 0),
                        c(0,  0,  0,  0,   0, 0, 0, 1),
                        c(0,  1,  0,  0,   0, 0, 0, 0),
                        c(0,  0,  1,  0,   0, 0, 0, 0))
        ci.opt <- c(-0.99999999, 1, 1, 0, 0)
      }

      if (weighting == "beta.unrestricted" & is.null(weighting.two) == FALSE) {
        if (weighting.two == "beta.restricted") {
          lf <- function(p) {
            llh_mf(df = df_llh,
                   y = ret,
                   x = covariate,
                   low.freq = low.freq,
                   mu = p["mu"],
                   omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                   alpha = p["alpha"], beta = p["beta"], gamma = p["gamma"],
                   m = p["m"], theta = p["theta"],
                   w1 = p["w1"], w2 = p["w2"], g_zero = g_zero, K = K,
                   x.two = covariate.two,
                   K.two = K.two, low.freq.two = low.freq.two,
                   theta.two = p["theta.two"], w1.two = 1, w2.two = p["w2.two"])
          }
          par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04,
                         m = 0, theta = 0, w1 = 1.00000001, w2 = 3, theta.two = 0, w2.two = 3)
          ui.opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0, 0, 0, 0),
                          c(0,  0,  0,    0, 0, 0, 1, 0, 0, 0),
                          c(0,  0,  0,    0, 0, 0, 0, 1, 0, 0),
                          c(0,  1,  0,    0, 0, 0, 0, 0, 0, 0),
                          c(0,  0,  1,    0, 0, 0, 0, 0, 0, 0),
                          c(0,  0,  0,    0, 0, 0, 0, 0, 0, 1))
          ci.opt <- c(-0.99999999, 1, 1, 0, 0, 1)
        }
      }
    }

    if (gamma == FALSE) {

      if (weighting == "beta.restricted") {
        lf <- function(p) {
          llh_mf(df = df_llh,
                 y = ret,
                 x = covariate,
                 low.freq = low.freq,
                 mu = p["mu"], omega = 1 - p["alpha"] - p["beta"],
                 alpha = p["alpha"], beta = p["beta"], gamma = 0,
                 m = p["m"], theta = p["theta"], w1 = 1, w2 = p["w2"],
                 g_zero = g_zero,
                 K = K)
        }
        par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, m = 0, theta = 0, w2 = 3)
        ui.opt <- rbind(c(0, -1, -1, 0, 0, 0),
                        c(0, 0, 0,  0, 0, 1),
                        c(0, 1, 0, 0,  0, 0),
                        c(0, 0, 1, 0, 0, 0))
        ci.opt <- c(-0.99999999, 1, 0, 0)
      }

      if (weighting == "beta.unrestricted") {
        lf <- function(p) {
          llh_mf(df = df_llh,
                 y = ret,
                 x = covariate,
                 low.freq = low.freq,
                 mu = p["mu"],
                 omega = 1 - p["alpha"] - p["beta"],
                 alpha = p["alpha"],
                 beta = p["beta"],
                 gamma = 0,
                 m = p["m"],
                 theta = p["theta"],
                 w1 = p["w1"],
                 w2 = p["w2"],
                 g_zero = g_zero,
                 K = K)
        }
        par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, m = 0, theta = 0, w1 = 1.00000001, w2 = 3)
        ui.opt <- rbind(c(0, -1, -1, 0, 0, 0, 0),
                        c(0,  0,  0, 0, 0, 1, 0),
                        c(0,  0,  0, 0, 0, 0, 1),
                        c(0,  1,  0, 0, 0, 0, 0),
                        c(0,  0,  1, 0, 0, 0, 0))
        ci.opt <- c(-0.99999999, 1, 1, 0, 0)
      }

    }

    if(is.null(control$par.start) == FALSE) {
      par.start <- control$par.start
    }

    p.e.nlminb <- constrOptim(theta = par.start, f = function(theta) { sum(lf(theta)) },
                              grad = NULL, ui = ui.opt, ci = ci.opt, hessian = FALSE)
    p.e.nlminb$value <- -p.e.nlminb$value

    if (multi.start == TRUE && gamma == TRUE) {
      p.e.nlminb.two <- try({
        suppressWarnings(optim(par = p.e.nlminb$par, fn = function (theta) {
          if( is.na(sum(lf(theta))) == TRUE | theta["alpha"] < 0 | theta["alpha"] + theta["beta"] + theta["gamma"]/2 >= 1 | theta["w2"] < 1) {
            NA
          } else {
            sum(lf(theta))
          }
        }, method = "BFGS"))}, silent = TRUE)

      if (class(p.e.nlminb.two) != "try-error" && -p.e.nlminb.two$value > p.e.nlminb$value) {
        p.e.nlminb <- p.e.nlminb.two
        p.e.nlminb$value <- -p.e.nlminb$value
      }

      par.max.lik.nr <- try({maxLik(logLik = function(x) - lf(x), start = par.start, method = "NR")}, silent = TRUE)
      if (class(par.max.lik.nr) != "try-error" && par.max.lik.nr$maximum > p.e.nlminb$value &&
          par.max.lik.nr$estimate["w2"] >= 1 &&
          par.max.lik.nr$estimate["alpha"] + par.max.lik.nr$estimate["beta"] + par.max.lik.nr$estimate["gamma"] / 2 < 1 &&
          par.max.lik.nr$estimate["alpha"] >= 0 && par.max.lik.nr$estimate["beta"] >= 0) {
        p.e.nlminb$par <- par.max.lik.nr$estimate
        p.e.nlminb$value <- par.max.lik.nr$maximum
      }
      par.max.lik.nm <- try({maxLik(logLik = function(x) -lf(x), start = par.start, method = "NM")}, silent = TRUE)
      if (class(par.max.lik.nm) != "try-error" && par.max.lik.nm$maximum > p.e.nlminb$value &&
          par.max.lik.nm$estimate["w2"] >= 1 &&
          par.max.lik.nm$estimate["alpha"] + par.max.lik.nm$estimate["beta"] + par.max.lik.nm$estimate["gamma"] / 2 < 1 &&
          par.max.lik.nm$estimate["alpha"] >= 0 && par.max.lik.nm$estimate["beta"] >= 0) {
        p.e.nlminb$par <- par.max.lik.nm$estimate
        p.e.nlminb$value <- par.max.lik.nm$maximum
      }

      p.e.nlminb.three <- try({
        suppressWarnings(optim(par = par.start, fn = function (theta) {
          if( is.na(sum(lf(theta))) == TRUE | theta["alpha"] < 0 | theta["alpha"] + theta["beta"] + theta["gamma"]/2  >= 1 | theta["w2"] < 1) {
            NA
          } else {
            sum(lf(theta))
          }
        }, method = "BFGS"))}, silent = TRUE)

      if (class(p.e.nlminb.three) != "try-error" && -p.e.nlminb.three$value > p.e.nlminb$value) {
        p.e.nlminb <- p.e.nlminb.three
        p.e.nlminb$value <- -p.e.nlminb$value
      }

    }

    if (multi.start == TRUE && gamma == FALSE) {

      p.e.nlminb.two <- try({
        suppressWarnings(optim(par = p.e.nlminb$par, fn = function (theta) {
          if( is.na(sum(lf(theta))) == TRUE | theta["alpha"] < 0 | theta["alpha"] + theta["beta"]  >= 1 | theta["w2"] < 1) {
            NA
          } else {
            sum(lf(theta))
          }
        }, method = "BFGS"))}, silent = TRUE)

      if (class(p.e.nlminb.two) != "try-error" && -p.e.nlminb.two$value > p.e.nlminb$value) {
        p.e.nlminb <- p.e.nlminb.two
        p.e.nlminb$value <- -p.e.nlminb$value
      }

      par.max.lik.nr <- try({maxLik(logLik = function(x) - lf(x), start = par.start, method = "NR")}, silent = TRUE)
      if (class(par.max.lik.nr) != "try-error" && par.max.lik.nr$maximum > p.e.nlminb$value &&
          par.max.lik.nr$estimate["w2"] >= 1 &&
          par.max.lik.nr$estimate["alpha"] + par.max.lik.nr$estimate["beta"] < 1 &&
          par.max.lik.nr$estimate["alpha"] >= 0 && par.max.lik.nr$estimate["beta"] >= 0) {
        p.e.nlminb$par <- par.max.lik.nr$estimate
        p.e.nlminb$value <- par.max.lik.nr$maximum
      }
      par.max.lik.nm <- try({maxLik(logLik = function(x) - lf(x), start = par.start, method = "NM")}, silent = TRUE)
      if (class(par.max.lik.nm) != "try-error" && par.max.lik.nm$maximum > p.e.nlminb$value &&
          par.max.lik.nm$estimate["w2"] >= 1 &&
          par.max.lik.nm$estimate["alpha"] + par.max.lik.nm$estimate["beta"] < 1 &&
          par.max.lik.nm$estimate["alpha"] >= 0 && par.max.lik.nm$estimate["beta"] >= 0) {
        p.e.nlminb$par <- par.max.lik.nm$estimate
        p.e.nlminb$value <- par.max.lik.nm$maximum
      }

      p.e.nlminb.three <- try({
        suppressWarnings(optim(par = par.start, fn = function (theta) {
          if( is.na(sum(lf(theta))) == TRUE | theta["alpha"] < 0 | theta["alpha"] + theta["beta"]  >= 1 | theta["w2"] < 1) {
            NA
          } else {
            sum(lf(theta))
          }
        }, method = "BFGS"))}, silent = TRUE)

      if (class(p.e.nlminb.three) != "try-error" && -p.e.nlminb.three$value > p.e.nlminb$value) {
        p.e.nlminb <- p.e.nlminb.three
        p.e.nlminb$value <- -p.e.nlminb$value
      }
    }
    par <- p.e.nlminb$par

    if (weighting == "beta.restricted") {
      if (is.null(x.two) == FALSE) {
        if (K.two > 1) {
          tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                                  w1 = 1, w2 = par["w2"], theta = par["theta"], m = par["m"], K = K,
                                  x.two = covariate.two, K.two = K.two, theta.two = par["theta.two"],
                                  low.freq.two = low.freq.two,
                                  w1.two = 1, w2.two = par["w2.two"])$tau
        } else {
          tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                                  w1 = 1, w2 = par["w2"], theta = par["theta"], m = par["m"], K = K,
                                  x.two = covariate.two, K.two = K.two, theta.two = par["theta.two"],
                                  low.freq.two = low.freq.two,
                                  w1.two = 1, w2.two = 1)$tau
        }

      } else {
        tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                                w1 = 1, w2 = par["w2"],
                                theta = par["theta"],
                                m = par["m"], K = K)$tau
      }

      tau_forecast <-
        exp(sum_tau_fcts(m = par["m"],
                         i = K + 1,
                         theta = par["theta"],
                         phivar = calculate_phi(w1 = 1, w2 = par["w2"], K = K),
                         covariate = c(tail(unlist(unique(data[c(x, low.freq)])[x]), K), NA),
                         K = K))

      if (is.null(x.two) == FALSE) {
        if (K.two > 1) {
          tau_forecast <-
            tau_forecast *
            exp(sum_tau_fcts(m = 0,
                             i = K.two + 1,
                             theta = par["theta.two"],
                             phivar = calculate_phi(w1 = 1, w2 = par["w2.two"], K = K.two),
                             covariate = c(tail(unlist(unique(data[c(x.two, low.freq.two)])[x.two]), K.two), NA),
                             K = K.two))
        } else {
          tau_forecast <-
            tau_forecast *
            exp(sum_tau_fcts(m = 0,
                             i = K.two + 1,
                             theta = par["theta.two"],
                             phivar = calculate_phi(w1 = 1, w2 = 1, K = K.two),
                             covariate = c(tail(unlist(unique(data[c(x.two, low.freq.two)])[x.two]), K.two), NA),
                             K = K.two))
        }

      }
    }

    if (weighting == "beta.unrestricted") {
      if (is.null(x.two) == FALSE) {
        tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                                w1 = par["w1"], w2 = par["w2"], theta = par["theta"], m = par["m"], K = K,
                                x.two = covariate.two, K.two = K.two, theta.two = par["theta.two"],
                                low.freq.two = low.freq.two,
                                w1.two = 1, w2.two = par["w2.two"])$tau
      } else {
        tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                                w1 = par["w1"], w2 = par["w2"],
                                theta = par["theta"],
                                m = par["m"], K = K)$tau
      }


      tau_forecast <-
        exp(sum_tau_fcts(m = par["m"],
                         i = K + 1,
                         theta = par["theta"],
                         phivar = calculate_phi(w1 = par["w1"], w2 = par["w2"], K = K),
                         covariate = c(tail(unlist(unique(data[c(x, low.freq)])[x]), K), NA),
                         K = K))
      if (is.null(x.two) == FALSE) {
        tau_forecast <-
          tau_forecast *
          exp(sum_tau_fcts(m = 0,
                           i = K.two + 1,
                           theta = par["theta.two"],
                           phivar = calculate_phi(w1 = 1, w2 = par["w2.two"], K = K.two),
                           covariate = c(tail(unlist(unique(data[c(x.two, low.freq.two)])[x.two]), K.two), NA),
                           K = K.two))
      }
    }

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

    if ((var.ratio.freq %in% c("date", low.freq)) == FALSE) {
      if (is.null(x.two) == TRUE) {
        df.fitted <- cbind(data[c("date", y, low.freq, x, var.ratio.freq)], g = g, tau = tau)
      } else {
        df.fitted <- cbind(data[c("date", y, low.freq, x, low.freq.two, x.two, var.ratio.freq)], g = g, tau = tau)
      }

    } else {
      if (is.null(x.two) == TRUE) {
        df.fitted <- cbind(data[c("date", y, low.freq, x)], g = g, tau = tau)
      } else {
        df.fitted <- cbind(data[c("date", y, low.freq, x, low.freq.two, x.two)], g = g, tau = tau)
      }

    }

    df.fitted$residuals <- unlist((df.fitted[y] - par["mu"]) / sqrt(df.fitted$g * df.fitted$tau))

  }
  df.fitted$date <- as.Date(date_backup)
  # Standard errors --------------------------------------------------------------------------------
  # inv_hessian <- try({
  #   solve(-optimHess(par = par, fn = function (theta) {
  #       if( is.na(sum(lf(theta))) == TRUE) {
  #         10000000
  #       } else {
  #         sum(lf(theta))
  #       }
  #     }))
  #   }, silent = TRUE)

  inv_hessian <- try({
    solve(-suppressWarnings(hessian(x = par, func = function (theta) {
        if( is.na(sum(lf(theta))) == TRUE) {
          0
        } else {
          -sum(lf(theta))
        }
      })))
    }, silent = TRUE)

  opg.std.err <- try({sqrt(diag(solve(crossprod(jacobian(func = function(theta) -lf(theta), x = par)))))},
                     silent = TRUE)
  if (class(opg.std.err)[1] == "try-error") {
    warning("Inverting the OPG matrix failed. No OPG standard errors calculated.")
    opg.std.err <- NA
  } else {
    opg.std.err <- opg.std.err * sqrt((mean(df.fitted$residuals^4, na.rm = TRUE) - 1) / 2)
  }

  if (class(inv_hessian)[1] == "try-error") {
    warning("Inverting the Hessian matrix failed. No robust standard errors calculated. Possible workaround: Multiply returns by 100.")
    rob.std.err <- NA
  } else {
    rob.std.err <- sqrt(diag(inv_hessian %*% crossprod(jacobian(func = lf, x = par)) %*% inv_hessian))
  }

  # Output -----------------------------------------------------------------------------------------
  output <-
    list(par = par,
         std.err = rob.std.err,
         broom.mgarch = data.frame(term = names(par),
                                   estimate = par,
                                   rob.std.err = rob.std.err,
                                   p.value = 2 * (1 - pnorm(unlist(abs(par/rob.std.err)))),
                                   opg.std.err = opg.std.err,
                                   opg.p.value = 2 * (1 - pnorm(unlist(abs(par/opg.std.err))))),
         tau = tau,
         g = g,
         df.fitted = df.fitted,
         K = K,
         weighting.scheme = weighting,
         llh = p.e.nlminb$value,
         bic = log(sum(!is.na(tau))) * length(par) - 2 * (p.e.nlminb$value),
         y = y,
         optim = p.e.nlminb)

  if (is.null(x.two) == FALSE) {
    output$K.two <- K.two
    output$weighting.scheme.two <- weighting.two
  }
  if (K == 0) {
    output$tau.forecast <- exp(par["m"])
  }


  # Additional output if there is a long-term component (K > 0) -------------------------------------
  if (K > 0) {
    output$variance.ratio <- 100 *
                     var(log(aggregate(df.fitted$tau, by = df.fitted[var.ratio.freq],
                                       FUN = mean)[,2]),
                         na.rm = TRUE) /
                     var(log(aggregate(df.fitted$tau * df.fitted$g, by = df.fitted[var.ratio.freq],
                                       FUN = mean)[,2]),
                         na.rm = TRUE)
    output$tau.forecast <- tau_forecast

    if (weighting == "beta.restricted") {
      output$est.weighting <- calculate_phi(1, w2 = par["w2"], K = K)
    }
    if (weighting == "beta.unrestricted") {
      output$est.weighting <- calculate_phi(w1 = par["w1"], w2 = par["w2"], K = K)
    }
    if (is.null(x.two) == FALSE) {
      if (K.two > 1) {
        output$est.weighting.two <- calculate_phi(w1 = 1, w2 = par["w2.two"], K = K.two)
      }
    }

  }

  # Add class mfGARCH for employing generic functions
  class(output) <- "mfGARCH"
  output
}
