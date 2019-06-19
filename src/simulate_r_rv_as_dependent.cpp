//Includes/namespaces
#include <Rcpp.h>
#include <cstdio>
using namespace Rcpp;

// [[Rcpp::export]]
List simulate_r_rv_as_dependent(double n_days, double n_intraday, double alpha, double beta, double gamma, NumericVector Z,
                double h0, int K, double m, double theta, NumericVector weights) {

  NumericVector r(n_days);
  NumericVector r_intraday(n_days * n_intraday);
  NumericVector h(n_days);
  NumericVector rv(n_days);
  NumericVector tau(n_days);
  NumericVector rv_tt(n_days);

  tau[0] = 1;
  rv_tt[0] = 1;
  h[0] = h0;

  for (int k = 0; k < n_intraday; k++) {
    r_intraday[k] = Z[k] * sqrt(h[0]) / sqrt(n_intraday);
  }
  r[0] = 0;
  rv[0] = 0;
  for (int k = 0; k < n_intraday; k++) {
    r[0] += r_intraday[k];
    rv[0] += (r_intraday[k] * r_intraday[k]);
  }


  for (int i = 1; i < n_days; i++) {
    // Rcout << "New day i: " << i << "\n";
    // Rcout << "The value of last-days return: " << r[i - 1] << "\n";
    if (r[i-1] < 0) {
      h[i] = 1 - alpha - beta - gamma / 2 + (alpha + gamma) * r[i-1] * r[i-1] / tau[i-1] + beta * h[i-1];
    } else {
      h[i] = 1 - alpha - beta - gamma / 2 + alpha * r[i-1] * r[i-1] / tau[i-1] + beta * h[i-1];
    }
    // Rcout << "Day i: " << i << "\n";
    // Rcout << "h: " << h[i] << "\n";
    rv_tt[i-1] = 0;
    for (int k = 0; k < 22; k++) {
      // Rcout << "The value of k: " << k << "\n";
      if (i - 1 - k >= 0) {
        // Rcout << "The value of returns_sq: " << r[i - i - k] * r[i - i - k] << "\n";
        // Rcout << "The value of rv22: " << rv_tt[i-1] << "\n";
        rv_tt[i-1] += r[i - 1 - k] * r[i - 1 - k] * 1 / 22;
        // Rcout << "The value of rv22: " << rv_tt[i-1] << "\n";
      }

    }
    // Rcout << "Day i: " << i << "\n";
    tau[i] = m;
    for (int j = 0; j < K; j++) {
      if (i - 1 - j >= 0) {
        tau[i] += theta * weights[j] * sqrt(rv_tt[i - 1 - j]);
      }
    }
    // Rcout << "Day i: " << i << "\n";
    tau[i] = exp(tau[i]);

    for (int k = (i * n_intraday); k < ((i+1) * n_intraday); k++) {
      r_intraday[k] = Z[k] * sqrt(h[i]) * sqrt(tau[i])/ sqrt(n_intraday);
    }
    // Rcout << "Day i: " << i << "\n";
    r[i] = 0;
    rv[i] = 0;
    for (int l = (i * n_intraday); l < ((i+1) * n_intraday); l++) {
      // Rcout << "Intraday l: " << l << "\n";
      r[i] += r_intraday[l];
      rv[i] += (r_intraday[l] * r_intraday[l]);
    }
  }

  return List::create(Named("ret_daily") = r,
                      Named("h_daily") = h,
                      Named("ret_intraday") = r_intraday,
                      Named("rv") = rv,
                      Named("rv22") = rv_tt,
                      Named("tau") = tau);
}
