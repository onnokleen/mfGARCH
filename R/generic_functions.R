#' @export
print.mfGARCH <- function(x) {
    if (class(x) != "mfGARCH") {
        stop("Obejct is not in class mfGARCH")
    } else {
        print("Print functions will be added soon.")
    }
}

# predict.mfGARCH <- function(x) {
#
# }

#' @export
predict.mfGARCH <- function(mgarch, horizon = c(1:10), fcts.tau = NULL, return = NULL, cond.var = NULL, cond.tau = NULL) {
  if (class(mgarch) != "mfGARCH") {
    stop("Obejct is not in class mfGARCH")
  }

  if (is.null(cond.var) == TRUE) {
    cond.var <- tail(mgarch$mgarch.g, 1)
  }

  if (is.null(cond.tau) == TRUE) {
    cond.tau <- tail(mgarch$mgarch.tau, 1)
  }

  if (is.null(fcts.tau) == TRUE) {
    fcts.tau <- mgarch$tau_forecast
  }

  if (is.null(return) == TRUE) {
    return <- tail(mgarch$df.fitted$return, 1)
  }

  fcts.tau * as.numeric(sapply(horizon, forecast_garch,
                               omega = 1 - mgarch$mgarch$par["alpha"] - mgarch$mgarch$par["beta"] - mgarch$mgarch$par["gamma"]/2,
                               alpha = mgarch$mgarch$par["alpha"],
                               beta = mgarch$mgarch$par["beta"],
                               gamma = mgarch$mgarch$par["gamma"],
                               ret = (return - - mgarch$par["mu"])/ sqrt(cond.tau),
                               g = cond.var))
}
