#' @export
print.mfGARCH <- function(x, ...) {
    if (class(x) != "mfGARCH") {
        stop("Obejct is not in class mfGARCH")
    } else {
        print(x$broom.mgarch)
    }
}

#' @importFrom graphics lines
#' @export
plot.mfGARCH <- function(x, ...) {
  if (class(x) != "mfGARCH") {
    stop("Obejct is not in class mfGARCH")
  }

  if (x$K == 0) {
    plot(x = x$df.fitted["date"], y = sqrt(x$df.fitted$g),
         type = "l",
         xlab = colnames(x$df.fitted)[3], ylab = "vol",
         main = "sqrt(g)", sub = "")
  } else {

    df_plot <- aggregate(x$df.fitted, by = list(x$df.fitted[, 3]), FUN = mean)
    plot(x = df_plot[, 1], y = sqrt(df_plot$g),
         type = "l",
         xlab = colnames(x$df.fitted)[3], ylab = "vol",
         main = "sqrt(tau * g) and sqrt(tau) in red", sub = "")
    lines(x = df_plot[, 1],
          y = sqrt(df_plot$tau),
          col = "red")
  }
}

#' @export
predict.mfGARCH <- function(object, horizon = c(1:10), fcts.tau = NULL, y.last = NULL, cond.var = NULL, cond.tau = NULL, ...) {
  if (class(object) != "mfGARCH") {
    stop("Obejct is not in class mfGARCH")
  }

  if (is.null(cond.var)) {
    cond.var <- tail(object$g, 1)
  }

  if (is.null(cond.tau)) {
    cond.tau <- tail(object$tau, 1)
  }

  if (is.null(fcts.tau)) {
    fcts.tau <- object$tau.forecast
  }

  if (is.null(y.last)) {
    return <- tail(object$df.fitted[object$y], 1)
  } else {
    return <- y.last
  }

  if (is.na(object$par["gamma"])) {
    fcts.tau * as.numeric(sapply(horizon, forecast_garch,
                                 omega = 1 - object$par["alpha"] - object$par["beta"] - 0/2,
                                 alpha = object$par["alpha"],
                                 beta = object$par["beta"],
                                 gamma = 0,
                                 ret = (return - object$par["mu"])/ sqrt(cond.tau),
                                 g = cond.var))
  } else {
    fcts.tau * as.numeric(sapply(horizon, forecast_garch,
                                 omega = 1 - object$par["alpha"] - object$par["beta"] - object$par["gamma"]/2,
                                 alpha = object$par["alpha"],
                                 beta = object$par["beta"],
                                 gamma = object$par["gamma"],
                                 ret = (return - object$par["mu"])/ sqrt(cond.tau),
                                 g = cond.var))
  }

}

#' This function plots the weighting scheme of an estimated GARCH-MIDAS model
#' @param x mfGARCH object obtained by fit_mfgarch
#' @importFrom graphics plot
#' @export
plot_weighting_scheme <- function(x) {
  if (class(x) != "mfGARCH") {
    stop("Obejct is not in class mfGARCH")
  }

  if (x$weighting.scheme == "beta.restricted") {
    k = c(1:x$K)
    phi = calculate_phi(w1 = 1, w2 = x$par["w2"], K = x$K)
    plot(k, phi)

  }

  if (x$weighting.scheme == "beta.unrestricted") {
    k = c(1:x$K)
    phi = calculate_phi(w1 = x$par["w1"], w2 = x$par["w2"], K = x$K)
    plot(k, phi)
  }

}
