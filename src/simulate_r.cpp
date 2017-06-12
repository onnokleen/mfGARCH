//Includes/namespaces
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List simulate_r(double n_days, double n_intraday, double alpha, double beta, double gamma, NumericVector Z, double h0) {

  NumericVector r(n_days);
  NumericVector r_intraday(n_days * n_intraday);
  NumericVector h(n_days);
  NumericVector rv(n_days);

  h[0] = h0;

  for (int k = 0; k < n_intraday; k++) {
    r_intraday[k] = Z[k] * sqrt(h[0]) / sqrt(n_intraday);
  }
  r[0] = 0;
  rv[0] = 0;
  for (int k = 0; k < n_intraday; k++) {
    r[0] = r[0] + r_intraday[k];
    rv[0] = rv[0] + (r_intraday[k] * r_intraday[k]);
  }


  for (int i = 1; i < n_days; i++) {
    if (r[i-1] < 0) {
      h[i] = 1 - alpha - beta - gamma / 2 + (alpha + gamma) * r[i-1] * r[i-1] + beta * h[i-1];
    } else {
      h[i] = 1 - alpha - beta - gamma / 2 + alpha * r[i-1] * r[i-1] + beta * h[i-1];
    }

    for (int k = (i * n_intraday); k < ((i+1) * n_intraday); k++) {
      r_intraday[k] = Z[k] * sqrt(h[i]) / sqrt(n_intraday);
    }

    r[i] = 0;
    rv[i] = 0;
    for (int l = (i * n_intraday); l < ((i+1) * n_intraday); l++) {
      r[i] = r[i] + r_intraday[l];
      rv[i] = rv[i] + (r_intraday[l] * r_intraday[l]);
    }
  }


  return List::create(Named("ret_daily") = r,
                      Named("h_daily") = h,
                      Named("ret_intraday") = r_intraday,
                      Named("rv") = rv);
}
