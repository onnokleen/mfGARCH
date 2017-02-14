#' This function estimates a multiplicative mixed-frequency GARCH model
#' @param df, y, x, K , low_freq, var_ratio_freq.
#' @keywords fit_mfgarch
#' @export fit_mfgarch
#' @importFrom magrittr %>%
#' @importFrom numDeriv jacobian
#' @importFrom stats nlminb
#' @importFrom stats optimHess
#' @importFrom tibble data_frame
# @examples likelihood_gjrgarch(0.01, 0.02, 0.9, 0.02, y = rnorm(1:4), mu = 0, g.0 = 0.2)

fit_mfgarch <- function(df, y, x, K, low_freq = "Date", var_ratio_freq = NULL) {

    if (sum(is.na(df$return) == TRUE) > 0 | sum(is.na(df$NFCI) == TRUE) > 0) {
        stop("Either y or x include NAs")
    }

    if ("Date" %in% colnames(df) == FALSE) {
        stop("No Date column.")
    }

    if (length(unlist(unique(select_(df, "Date")))) != dim(df)[1]) {
        stop("There is more than one observation per high frequency (presumably day).")
    }

    # If this is not the case, we may order by low frequency variable
    df <- df %>% dplyr::arrange_("Date")

    if (y %in% colnames(df) == FALSE) {
        stop(paste("There is no variable in your data frame with name ", y, "."))
    }
    if (x %in% colnames(df) == FALSE) {
        stop(paste("There is no variable in your data frame with name ", x, "."))
    }
    if (low_freq %in% colnames(df) == FALSE) {
        stop(paste("There is no low freq. variable in your data frame with name ", low_freq, "."))
    }

    if (is.null(var_ratio_freq) == TRUE) {
        var_ratio_freq <- low_freq
        print("No frequency specified for calculating the variance ratio - default: low_freq")
    }

    df.llh <- df %>% select_(., ~get(y), ~get(x), ~get(low_freq))

    g_zero <- var(df$return^2)

    # Parameter estimation
    if (K == 0) {
        lf <- function(parameters) {
            sum(likelihood_mg_simple(df = df.llh,
                                     y = y,
                                     mu = parameters["mu"],
                                     omega = 1 - parameters["alpha"] - parameters["beta"] - parameters["gamma"]/2,
                                     alpha = parameters["alpha"],
                                     beta = parameters["beta"],
                                     gamma = parameters["gamma"],
                                     m = parameters["m"],
                                     g_zero = g_zero))
        }


        parameters.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04, m = 0)

        ui.opt <- rbind(c(0, -1, -1, -1/2, 0), c(0, 1, 0, 0, 0), c(0, 0, 1, 0, 0))
        ci.opt <- c(-0.99999, 0, 0)


        p.e.nlminb <- constrOptim(theta = parameters.start, f = lf, grad = NULL, ui = ui.opt, ci = ci.opt, hessian = FALSE)
        returns <- df %>% select_(~get(y)) %>% unlist() %>% as.numeric()
        tau.estimate <- rep(exp(p.e.nlminb$par["m"]), times = length(returns))
        g.estimate.mg <- c(rep(NA, times = sum(is.na((returns - p.e.nlminb$par["mu"])/sqrt(tau.estimate)))), calculate_g(omega = 1 - p.e.nlminb$par["alpha"] - p.e.nlminb$par["beta"] -
            p.e.nlminb$par["gamma"]/2, alpha = p.e.nlminb$par["alpha"], beta = p.e.nlminb$par["beta"], gamma = p.e.nlminb$par["gamma"], as.numeric(na.exclude((returns - p.e.nlminb$par["mu"])/sqrt(tau.estimate))),
            g0 = g_zero))

        df.fitted <- data_frame(returns = returns, g.mgarch = g.estimate.mg, tau.mgarch = tau.estimate) %>% mutate(., residuals = (returns - p.e.nlminb$par["mu"])/sqrt(g.mgarch *
            tau.mgarch))
    }

    if (K == 1) {
        lf <- function(parameters) {
            sum(likelihood_mg_mf(df = df.llh, y = y, x = x, low_freq = low_freq, mu = parameters["mu"], omega = 1 - parameters["alpha"] - parameters["beta"] - parameters["gamma"]/2, alpha = parameters["alpha"],
                beta = parameters["beta"], gamma = parameters["gamma"], m = parameters["m"], theta = parameters["theta"], w1 = 1, w2 = 1, g_zero = g_zero, K = K))
        }

        parameters.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04, m = 0, theta = 0)

        ui.opt <- rbind(c(0, -1, -1, -1/2, 0, 0), c(0, 1, 0, 0, 0, 0), c(0, 0, 1, 0, 0, 0))
        ci.opt <- c(-0.99999, 0, 0)

        p.e.nlminb <- constrOptim(theta = parameters.start, f = lf, grad = NULL, ui = ui.opt, ci = ci.opt, hessian = FALSE)

        tau.estimate <- calculate_tau_mf(df = df, x = x, low_freq = low_freq, w1 = 1, w2 = 1, theta = p.e.nlminb$par["theta"], m = p.e.nlminb$par["m"], K = K)$tau

        returns <- df$return

        g.estimate.mg <- c(rep(NA, times = sum(is.na((returns - p.e.nlminb$par["mu"])/sqrt(tau.estimate)))), calculate_g(omega = 1 - p.e.nlminb$par["alpha"] - p.e.nlminb$par["beta"] -
            p.e.nlminb$par["gamma"]/2, alpha = p.e.nlminb$par["alpha"], beta = p.e.nlminb$par["beta"], gamma = p.e.nlminb$par["gamma"], as.numeric(na.exclude((returns - p.e.nlminb$par["mu"])/sqrt(tau.estimate))),
            g0 = g_zero))

        df.fitted <- data_frame(returns = returns, g.mgarch = g.estimate.mg, tau.mgarch = tau.estimate) %>% mutate(., residuals = (returns - p.e.nlminb$par["mu"])/sqrt(g.mgarch *
            tau.mgarch))

    }

    if (K > 1) {
        lf <- function(parameters) {
            sum(likelihood_mg_mf(df = df.llh, y = y, x = x, low_freq = low_freq, mu = parameters["mu"], omega = 1 - parameters["alpha"] - parameters["beta"] - parameters["gamma"]/2,
                alpha = parameters["alpha"], beta = parameters["beta"], gamma = parameters["gamma"], m = parameters["m"], theta = parameters["theta"], w1 = 1, w2 = parameters["w2"],
                g_zero = g_zero, K = K))
        }

        parameters.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04, m = 0, theta = 0, w2 = 3)

        ui.opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0), c(0, 0, 0, 0, 0, 0, 1), c(0, 1, 0, 0, 0, 0, 0), c(0, 0, 1, 0, 0, 0, 0))
        ci.opt <- c(-0.99999999, 1, 0, 0)

        p.e.nlminb <- constrOptim(theta = parameters.start, f = lf, grad = NULL, ui = ui.opt, ci = ci.opt, hessian = FALSE)

        tau.estimate <- calculate_tau_mf(df = df, x = x, low_freq = low_freq, w1 = 1, w2 = p.e.nlminb$par["w2"], theta = p.e.nlminb$par["theta"], m = p.e.nlminb$par["m"], K = K)$tau

        returns <- df$return

        g.estimate.mg <- c(rep(NA, times = sum(is.na((returns - p.e.nlminb$par["mu"])/sqrt(tau.estimate)))), calculate_g(omega = 1 - p.e.nlminb$par["alpha"] - p.e.nlminb$par["beta"] -
            p.e.nlminb$par["gamma"]/2, alpha = p.e.nlminb$par["alpha"], beta = p.e.nlminb$par["beta"], gamma = p.e.nlminb$par["gamma"], as.numeric(na.exclude((returns - p.e.nlminb$par["mu"])/sqrt(tau.estimate))),
            g0 = g_zero))

        df.fitted <- df %>% cbind(g.mgarch = g.estimate.mg, tau.mgarch = tau.estimate) %>% mutate(., residuals = (returns - p.e.nlminb$par["mu"])/sqrt(g.mgarch * tau.mgarch))
        # data_frame(returns = returns, g.mgarch = g.estimate.mg, tau.mgarch = tau.estimate) %>%

    }

    p.e.nlminb.n.i.hess <- solve(-optimHess(par = p.e.nlminb$par, fn = lf))

    if (K == 0) {
        lf_mgarch_score <- function(parameters) {
            likelihood_mg_simple(df = df.llh, y = y, mu = parameters["mu"], omega = 1 - parameters["alpha"] - parameters["beta"] - parameters["gamma"]/2, alpha = parameters["alpha"],
                beta = parameters["beta"], gamma = parameters["gamma"], m = parameters["m"], g_zero = g_zero)
        }
    }

    if (K == 1) {
        lf_mgarch_score <- function(parameters) {
            likelihood_mg_mf(df = df.llh, y = y, x = x, low_freq = low_freq, mu = parameters["mu"], omega = 1 - parameters["alpha"] - parameters["beta"] - parameters["gamma"]/2, alpha = parameters["alpha"],
                beta = parameters["beta"], gamma = parameters["gamma"], m = parameters["m"], theta = parameters["theta"], w1 = 1, w2 = 1, g_zero = g_zero, K = K)
        }
    }

    if (K > 1) {
        lf_mgarch_score <- function(parameters) {
            likelihood_mg_mf(df = df.llh, y = y, x = x, low_freq = low_freq, mu = parameters["mu"], omega = 1 - parameters["alpha"] - parameters["beta"] - parameters["gamma"]/2, alpha = parameters["alpha"],
                beta = parameters["beta"], gamma = parameters["gamma"], m = parameters["m"], theta = parameters["theta"], w1 = 1, w2 = parameters["w2"], g_zero = g_zero, K = K)
        }
    }

    G <- numDeriv::jacobian(func = lf_mgarch_score, x = p.e.nlminb$par)

    GGsum <- t(G) %*% G
    # p.e.nlminb.robust.standard.errors <- sqrt(diag(solve(GGsum))) * mean((df.fitted$residuals^2 - 1)^2, na.rm = TRUE)/2

    p.e.nlminb.robust.standard.errors <- (p.e.nlminb.n.i.hess %*% GGsum %*% p.e.nlminb.n.i.hess) %>% diag() %>% sqrt()

    if (K == 0) {
        return(list(mgarch = p.e.nlminb, broom.mgarch = data_frame(term = c("mu", "alpha", "beta", "gamma", "m", "llh"), estimate = c(p.e.nlminb$par[1:5], -p.e.nlminb$value), rob.std.err = c(p.e.nlminb.robust.standard.errors[1:5],
            NA)), mgarch.tau = tau.estimate, mgarch.g = g.estimate.mg, df.fitted = df.fitted))
    }

    if (K == 1) {
        return(list(mgarch = p.e.nlminb, broom.mgarch = data_frame(term = c("mu", "alpha", "beta", "gamma", "m", "theta", "llh"), estimate = c(p.e.nlminb$par[1:6], -p.e.nlminb$value),
            rob.std.err = c(p.e.nlminb.robust.standard.errors[1:6], NA)), mgarch.tau = tau.estimate, mgarch.g = g.estimate.mg, df.fitted = df.fitted))

    }

    if (K > 1) {
        return(list(mgarch = p.e.nlminb, broom.mgarch = data_frame(term = c("mu", "alpha", "beta", "gamma", "m", "theta", "w2", "llh"), estimate = c(p.e.nlminb$par[1:7], -p.e.nlminb$value),
            rob.std.err = c(p.e.nlminb.robust.standard.errors[1:7], NA)), mgarch.tau = tau.estimate, mgarch.g = g.estimate.mg, df.fitted = df.fitted, variance.ratio = 100 * (df.fitted %>%
            group_by_(var_ratio_freq) %>% summarise(mean.tau = mean(tau.mgarch), mean.tau.g = mean(tau.mgarch * g.mgarch)) %>% ungroup() %>% summarise(var.ratio = var(log(mean.tau),
            na.rm = TRUE)/var(log(mean.tau.g), na.rm = TRUE)) %>% unlist())))
    }
}


calculate_tau_mf <- function(df, x, low_freq, w1, w2, theta, m, K) {

    phi.var <- calculate_phi(w1, w2, K)

    covariate <- df %>% select_(., ~get(low_freq), ~get(x)) %>% unique() %>% select_(., ~get(x)) %>% unlist()

    covariate <- c(rep(NA, times = K), covariate)
    tau <- exp(sapply(c((K + 1):length(covariate)), FUN = sum_tau, m = m, theta = theta, phivar = phi.var, covariate = covariate, K = K))
    full_join(df, eval(parse(text = paste0("tbl_df(cbind(", low_freq, " = unlist(unique(select_(df, ~get(low_freq)))), tau))"))), by = low_freq)

}

likelihood_mg_mf <- function(df, x, y, low_freq, mu, omega, alpha, beta, gamma, m, theta, w1 = 1, w2 = 1, g_zero, K = 2) {

    tau <- calculate_tau_mf(df = df, x = x, low_freq = low_freq, w1 = w1, w2 = w2, theta = theta, m = m, K = K)$tau

    ret <- df %>% select_(., ~get(y)) %>% unlist()

    ret <- ret[which.min(is.na(tau)):length(ret)]  # lags can't be used for likelihood
    tau <- tau[which.min(is.na(tau)):length(tau)]

    g <- calculate_g(omega = omega, alpha = alpha, beta = beta, gamma = gamma, returns = ((ret - mu)/sqrt(tau)), g0 = g_zero)
    if (sum(g <= 0) > 0) {
        rep(NA, times = length(y))
    } else {
        1/2 * log(2 * pi) + 1/2 * log(g * tau) + 1/2 * (ret - mu)^2/(g * tau)
    }
}


likelihood_mg_simple <- function(df, y, mu, omega, alpha, beta, gamma, m, g_zero) {

    ret <- df %>% select_(., ~get(y)) %>% unlist()

    tau <- rep(exp(m), times = length(ret))

    g <- calculate_g(omega = omega, alpha = alpha, beta = beta, gamma = gamma, returns = ((ret - mu)/sqrt(exp(m))), g0 = g_zero)

    1/2 * log(2 * pi) + 1/2 * log(g * tau) + 1/2 * (ret - mu)^2/(g * tau)
}





