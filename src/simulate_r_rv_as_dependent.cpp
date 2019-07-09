//Includes/namespaces
#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List simulate_r_rv_as_dependent(double n_days, double n_intraday, double alpha, double beta, double gamma, NumericVector Z,
                                  double h0, int K, double m, double theta, NumericVector weights, int lowfreq, bool rvol) {

  NumericVector r(n_days);
  NumericVector r_intraday(n_days * n_intraday);
  NumericVector h(n_days);
  NumericVector rv(n_days);
  NumericVector taulowfreq(n_days/lowfreq);
  NumericVector rv_tt(n_days/lowfreq);

  taulowfreq[0] = exp(m);
  rv_tt[0] = 1;
  h[0] = h0;

  for (int k = 0; k < n_intraday; k++) {
    r_intraday[k] = Z[k] * std::sqrt(h[0]) / std::sqrt(n_intraday);
  }
  r[0] = 0;
  rv[0] = 0;
  for (int k = 0; k < n_intraday; k++) {
    r[0] += r_intraday[k];
    rv[0] += (r_intraday[k] * r_intraday[k]);
  }


  for (int i = 1; i < n_days; i++) {
    // Rcout << "New day i: " << i << "\n" << std::endl;
    // Rcout << "The value of last-days return: " << r[i - 1] << "\n";
    if (r[i-1] < 0) {
      h[i] = 1 - alpha - beta - gamma / 2 + (alpha + gamma) * pow(r[i-1], 2) / taulowfreq[(i -1)/lowfreq] + beta * h[i-1];
    } else {
      h[i] = 1 - alpha - beta - gamma / 2 + alpha * pow(r[i-1] / r[i-1], 2) / taulowfreq[(i - 1)/lowfreq] + beta * h[i-1];
    }
    // Rcpp::Rcout << "Day i: " << i << "\n";
    // Rcout << "h: " << h[i] << "\n";
    // rv_tt[i-1] = 0;
    // for (int k = 0; k < 22; k++) {
    //   // Rcout << "The value of k: " << k << "\n";
    //   if (i - 1 - k >= 0) {
    //     // Rcout << "The value of returns_sq: " << r[i - i - k] * r[i - i - k] << "\n";
    //     // Rcout << "The value of rv22: " << rv_tt[i-1] << "\n";
    //     rv_tt[i-1] += r[i - 1 - k] * r[i - 1 - k] * 1 / 22;
    //     // Rcout << "The value of rv22: " << rv_tt[i-1] << "\n";
    //   }
    //
    // }
    // Rcout << "Day i: " << i << "\n";
    if (i%lowfreq == 0 && i/lowfreq <= n_days/lowfreq) {
      rv_tt[i/lowfreq - 1] = 0;
      for (int j = 1; j <= lowfreq; j++) {
        if (i - j >= 0) {
          rv_tt[i / lowfreq - 1] += pow(r[i - j], 2);

        }
      }
      // Rcout << "rv_tt: " << rv_tt[(i+1)/lowfreq - 1] << "\n";
      if (rvol == TRUE) {
        rv_tt[i/lowfreq - 1] = std::sqrt(rv_tt[i/lowfreq - 1] / lowfreq);
      }

    }

    if (i%lowfreq == 0 && i/lowfreq <= n_days/lowfreq) {
      taulowfreq[i/lowfreq] = m;
      for (int j = 0; j < K; j++) {
        if (i/lowfreq - 1 - j >= 0) {
          taulowfreq[i/lowfreq] += theta * weights[j] * rv_tt[i/lowfreq - 1 - j];
          // Rcout << "rv_tt in taulowfreq: " << sqrt(rv_tt[(i+1)/lowfreq - 1 - j]) << "\n";
        }
      }
      taulowfreq[i/lowfreq] = std::exp(taulowfreq[i/lowfreq]);
    }



    for (int k = (i * n_intraday); k < ((i+1) * n_intraday); k++) {
      r_intraday[k] = Z[k] * std::sqrt(h[i] * taulowfreq[i/lowfreq] / n_intraday);
    }
    // Rcout << "Day i: " << i << "\n";
    r[i] = 0;
    rv[i] = 0;
    for (int l = (i * n_intraday); l < ((i+1) * n_intraday); l++) {
      // Rcout << "Intraday l: " << l << "\n";
      r[i] += r_intraday[l];
      rv[i] += pow(r_intraday[l], 2);
    }
  }

  // rv_tt = std::sqrt(rv_tt / 22);

  return List::create(Named("ret_daily") = r,
                      Named("h_daily") = h,
                      Named("ret_intraday") = r_intraday,
                      Named("rv") = rv,
                      Named("rv22") = rv_tt,
                      Named("tau") = taulowfreq);
}
