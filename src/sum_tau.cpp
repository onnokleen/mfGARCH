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
//' @export
// [[Rcpp::export]]

double sum_tau(int i, double m, double theta, NumericVector phivar, NumericVector covariate, int K) {

  double exponential = m;

  for (int j = 1; j <= K; j++) {
    exponential += theta * phivar[j-1] * covariate[i - 1 - j];
  }

  return exponential;
}





