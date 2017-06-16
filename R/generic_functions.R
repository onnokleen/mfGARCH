#' @export
print.mfGARCH <- function(x, ...) {
    if (class(x) != "mfGARCH") {
        stop("Obejct is not in class mfGARCH")
    } else {
        print(x$broom.mgarch)
    }
}

# predict.mfGARCH <- function(x) {
#
# }

#' @export
predict.mfGARCH <- function(object, horizon = c(1:10), fcts.tau = NULL, return = NULL, cond.var = NULL, cond.tau = NULL, ...) {
  if (class(object) != "mfGARCH") {
    stop("Obejct is not in class mfGARCH")
  }

  if (is.null(cond.var) == TRUE) {
    cond.var <- tail(object$mgarch.g, 1)
  }

  if (is.null(cond.tau) == TRUE) {
    cond.tau <- tail(object$mgarch.tau, 1)
  }

  if (is.null(fcts.tau) == TRUE) {
    fcts.tau <- object$tau_forecast
  }

  if (is.null(return) == TRUE) {
    return <- tail(object$df.fitted$return, 1)
  }

  fcts.tau * as.numeric(sapply(horizon, forecast_garch,
                               omega = 1 - object$mgarch$par["alpha"] - object$mgarch$par["beta"] - object$mgarch$par["gamma"]/2,
                               alpha = object$mgarch$par["alpha"],
                               beta = object$mgarch$par["beta"],
                               gamma = object$mgarch$par["gamma"],
                               ret = (return - - object$par["mu"])/ sqrt(cond.tau),
                               g = cond.var))
}
