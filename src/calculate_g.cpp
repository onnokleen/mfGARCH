
#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
//' @export calculate_g
// [[Rcpp::export]]
NumericVector calculate_g(double omega, double alpha, double beta, double gamma, NumericVector returns, double g0) {
  int n = returns.size();

  NumericVector g(n);
  g[0] = g0;

//  NumericVector constant(n);
//  constant[0] = omega;

  for (int i = 1; i < n; i++) {
//    constant[i] = omega;

    if (returns[i-1] >= 0) {
      g[i] = omega + alpha * returns[i-1] * returns[i-1] + beta * g[i-1];
    } else {
      g[i] = omega + alpha * returns[i-1] * returns[i-1] + gamma * returns[i-1] * returns[i-1] + beta * g[i-1];
    }
  }

  return g;
}

