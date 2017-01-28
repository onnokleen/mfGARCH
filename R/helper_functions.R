# Calculates the phi weighting scheme
calculate_phi <- function(w1, w2, K) {
  weights <- sapply(c(1:K),
                    FUN = function(j) (j / (K + 1)) ^ (w1 - 1) * (1 - j / (K + 1)) ^ (w2 - 1))
  weights <- weights / sum(weights)
  weights
}

# Calculates the long-term component in its log specification
calculate_tau <- function(covariate, w1, w2, theta, m, K) {
  phi_var <- calculate_phi(w1, w2, K)
  covariate <- c(rep(NA, times = K), covariate)
  ln_tau <- exp(
    sapply(c( (K + 1):length(covariate) ),
           FUN = sum_tau,
           m = m, theta = theta, phivar = phi_var, covariate = covariate, K = K)
    )
  ln_tau
}
