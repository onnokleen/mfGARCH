#' @keywords internal
forecast_garch <- function(omega, alpha, beta, gamma, g, ret, steps.ahead) {
  omega / (1 - alpha - gamma/2 - beta) + (alpha + beta + gamma/2)^(steps.ahead - 1) * (omega + (alpha + gamma/2 * as.numeric(ret < 0)) * ret^2 + beta * g - omega / (1 - alpha - gamma/2 - beta))
}
