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
    cond.var <- tail(object$g, 1)
  }

  if (is.null(cond.tau) == TRUE) {
    cond.tau <- tail(object$tau, 1)
  }

  if (is.null(fcts.tau) == TRUE) {
    fcts.tau <- object$tau_forecast
  }

  if (is.null(return) == TRUE) {
    return <- tail(object$df.fitted$return, 1)
  }

  if (is.na(object$par["gamma"]) == TRUE) {
    fcts.tau * as.numeric(sapply(horizon, forecast_garch,
                                 omega = 1 - object$par["alpha"] - object$par["beta"] - 0/2,
                                 alpha = object$par["alpha"],
                                 beta = object$par["beta"],
                                 gamma = 0,
                                 ret = (return - - object$par["mu"])/ sqrt(cond.tau),
                                 g = cond.var))
  } else {
    fcts.tau * as.numeric(sapply(horizon, forecast_garch,
                                 omega = 1 - object$par["alpha"] - object$par["beta"] - object$par["gamma"]/2,
                                 alpha = object$par["alpha"],
                                 beta = object$par["beta"],
                                 gamma = object$par["gamma"],
                                 ret = (return - - object$par["mu"])/ sqrt(cond.tau),
                                 g = cond.var))
  }

}

#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 ggplot
plot_weighting_scheme <- function(x) {
  if (class(x) != "mfGARCH") {
    stop("Obejct is not in class mfGARCH")
  }

  if (x$weighting.scheme == "beta.one.sided") {
    df_weighting <-
      data_frame(
        k = c(1:x$K),
        phi = calculate_phi(w1 = 1, w2 = x$par["w2"], K = x$K))
  }

  if (x$weighting.scheme == "beta.two.sided") {
    df_weighting <-
      data_frame(k = c(1:x$K),
                 phi = calculate_phi(w1 = x$par["w1"], w2 = x$par["w2"], K = x$K))
  }

  ggplot() +
    geom_line(data = df_weighting, aes(x = k, y = phi))

}
