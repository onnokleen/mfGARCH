//Includes/namespaces
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector calculate_g(double omega, double alpha, double beta, double gamma, NumericVector returns, double g0) {
  int n = returns.size();

  NumericVector g(n);
  g[0] = g0;

  for (int i = 1; i < n; i++) {

    if (returns[i-1] >= 0) {
      g[i] = omega + alpha * returns[i-1] * returns[i-1] + beta * g[i-1];
    } else {
      g[i] = omega + alpha * returns[i-1] * returns[i-1] + gamma * returns[i-1] * returns[i-1] + beta * g[i-1];
    }
  }

  return g;
}
