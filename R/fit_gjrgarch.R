#' @importFrom magrittr %>%
#' @importFrom tibble data_frame

fit_gjrgarch <- function(y) {

    lf_garch <- function(parameters) {
        sum(likelihood_gjrgarch(mu = parameters["mu"], omega = parameters["omega"], alpha = parameters["alpha"], beta = parameters["beta"], gamma = parameters["gamma"], y = returns,
            g.0 = parameters["g.0"]), na.rm = TRUE)
    }



    p.e.start <- c(mu = 0, omega = 0.02, alpha = 0.04, beta = 0.85, gamma = 0, g.0 = 0.2)

    # Estimate parameters
    p.e <- stats::nlminb(start = p.e.start, objective = lf_garch, lower = c(-Inf, 0.001, 0.001, 0.001, -Inf, 1e-04), control = list(iter.max = 500, eval.max = 400))

    # Calculate Hessian
    p.e.hess <- stats::optimHess(par = p.e$par, fn = lf_garch)

    # Score function - derivative w.r.t. g.0 is not interesting
    lf_garch_score <- function(parameters) {
        likelihood_gjrgarch(mu = parameters["mu"], omega = parameters["omega"], alpha = parameters["alpha"], beta = parameters["beta"], gamma = 0, y = returns, g.0 = p.e$par["g.0"])
    }


    jacobian.score <- numDeriv::jacobian(func = lf_garch_score, x = p.e$par[1:4])

    sum.of.out.prod.score <- t(jacobian.score) %*% jacobian.score

    p.e.gjr.robust.standard.errors <- (solve(-p.e.hess[1:4, 1:4]) %*% sum.of.out.prod.score %*% solve(-p.e.hess[1:4, 1:4])) %>% diag() %>% sqrt()


    # Calculate fitted values

    g.estimate <- grch::calculate_g(omega = p.e$par["omega"], alpha = p.e$par["alpha"], beta = p.e$par["beta"], p.e$par["gamma"], returns - p.e$par["mu"], p.e$par["g.0"])
    augment <- data_frame(returns = returns, g.gjr = g.estimate)
    # Data frame for parameters
    tidy <- data_frame(term = c("mu", "omega", "alpha", "beta"), estimate = p.e$par[1:4], rob.std.err = p.e.gjr.robust.standard.errors)
    # Data frame for objectives of grch-object
    glance <- data_frame(LLH = p.e$objective)

    list(tidy, augment, glance)
}
