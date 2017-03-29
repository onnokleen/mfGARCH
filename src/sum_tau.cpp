#include <Rcpp.h>
using namespace Rcpp;

//' @export sum_tau
// [[Rcpp::export(name = "sum_tau")]]
// Calculates the long-term component
double sum_tau(int i, double m, double theta, NumericVector phivar, NumericVector covariate, int K) {

  double exponential = m;

  for (int j = 1; j <= K; j++) {
    exponential += theta * phivar[j-1] * covariate[i - 1 - j];
  }

  return exponential;
}
