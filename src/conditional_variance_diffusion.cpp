#include <Rcpp.h>
using namespace Rcpp;
// 
// // [[Rcpp::export]]
// List simulate_r(double n_days, double n_intraday, double alpha, double beta, double gamma, NumericVector Z, double h0) {
//   int n = n_days * n_intraday;
//   NumericVector r(n_days);
//   NumericVector r_intraday(n_days * n_intraday);
//   NumericVector h(n_days);
//   NumericVector rv(n_days);
//   
//   h[0] = h0;
//   
//   for (int k = 0; k < n_intraday; k++) {
//     r_intraday[k] = Z[k] * sqrt(h[0]) / sqrt(n_intraday);
//   }
//   r[0] = 0;
//   rv[0] = 0;
//   for (int k = 0; k < n_intraday; k++) {
//     r[0] = r[0] + r_intraday[k];
//     rv[0] = rv[0] + (r_intraday[k] * r_intraday[k]);
//   }
// 
//   
//   for (int i = 1; i < n_days; i++) {
//     if (r[i-1] < 0) {
//       h[i] = 1 - alpha - beta - gamma / 2 + (alpha + gamma) * r[i-1] * r[i-1] + beta * h[i-1];
//     } else {
//       h[i] = 1 - alpha - beta - gamma / 2 + alpha * r[i-1] * r[i-1] + beta * h[i-1];
//     }
// 
//     for (int k = (i * n_intraday); k < ((i+1) * n_intraday); k++) {
//       r_intraday[k] = Z[k] * sqrt(h[i]) / sqrt(n_intraday);
//     }
//     
//     r[i] = 0;
//     rv[i] = 0;
//     for (int l = (i * n_intraday); l < ((i+1) * n_intraday); l++) {
//       r[i] = r[i] + r_intraday[l];
//       rv[i] = rv[i] + (r_intraday[l] * r_intraday[l]);
//     }
//   }
//   
//   
//   return List::create(Named("ret_daily") = r,
//                       Named("h_daily") = h,
//                       Named("ret_intraday") = r_intraday,
//                       Named("rv") = rv);
// }

// // [[Rcpp::export]]
// NumericVector calculate_h(double ndays, double delta, double mu, double alpha, double beta, double gamma, NumericVector Z1, NumericVector Z2, double pi, double h0) {
//   int n = ndays * delta;
//   NumericVector h(n);
//   h[0] = h0;
//   for (int i = 1; i < n; i++) {
//     h[i] = h[i-1] + (1 - alpha - beta - gamma/2 + (beta + alpha + gamma/2 - 1) * h[i-1]) * 1 / delta - gamma * sqrt(2/pi) * h[i-1] * sqrt(1 / delta) * Z1[i-1] + sqrt(2 * alpha * alpha + (1 * pi - 2) / (1 * pi) * gamma * gamma + 2 * alpha * gamma) * h[i-1] * sqrt(1 / delta) * Z2[i-1];
//   }
//   return h;
// }

// // [[Rcpp::export]]
// NumericVector calculate_h(double ndays, double delta, double mu, double alpha, double beta, double gamma, NumericVector Z1, NumericVector Z2, double pi, double h0) {
//   int n = ndays * delta;
//   NumericVector h(n);
//   h[0] = h0;
//   for (int i = 1; i < n; i++) {
//     h[i] = h[i-1] + (1 - alpha - beta - gamma/2 + (beta + alpha + gamma/2 - 1) * h[i-1]) * 1 / delta - gamma * sqrt(2/pi) * h[i-1] * sqrt(1 / delta) * Z1[i-1] + sqrt(2 * alpha * alpha + (1 * pi - 2) / (1 * pi) * gamma * gamma + 2 * alpha * gamma) * h[i-1] * sqrt(1 / delta) * Z2[i-1];
//   }
//   return h;
// }
// 
// 
// [[Rcpp::export]]
NumericVector calculate_h_andersen(double ndays, double delta, double mu, double theta, double omega, double lambda, NumericVector Z, double pi, double h0) {
  int n = ndays * delta;
  NumericVector h(n);
  h[0] = h0;
  for (int i = 1; i < n; i++) {
    //h[i] = h[i-1] + theta * omega * 1 / delta - h[i-1] * theta * 1/delta + sqrt(2 * lambda * theta / delta) * h[i-1] * Z1[i-1];
    h[i] = theta * omega * 1/delta + h[i-1] * (1 - theta * 1/delta + sqrt(2 * lambda * theta * 1/delta) * Z[i-1]);
  }
  return h;
}
// 
// NumericVector calculate_h_gjr(double ndays, double delta, double theta, double omega, double lambda, double gamma, NumericVector Z,  double h0) {
//   int n = ndays * delta;
//   NumericVector h(n);
//   h[0] = h0;
//   for (int i = 1; i < n; i++) {
//     //h[i] = h[i-1] + theta * omega * 1 / delta - h[i-1] * theta * 1/delta + sqrt(2 * lambda * theta / delta) * h[i-1] * Z1[i-1];
//     h[i] = theta * omega * 1/delta + h[i-1] * (1 - theta * 1/delta + sqrt(2 * lambda * theta + 0.75 * gamma * gamma) * sqrt(1/delta) * Z[i-1]);
//   }
//   return h;
// }
// 
// [[Rcpp::export]]
NumericVector calculate_p(double ndays, double delta, double mu, NumericVector Zp, NumericVector h, double p0) {
  int n = ndays * delta;
  NumericVector p(n);
  p[0] = p0;
  for (int i = 1; i < n; i++) {
    p[i] = p[i-1] + sqrt(h[i]) * sqrt(1 / delta) * Zp[i] ;
    //h[i] = theta * omega * 1 / delta + h[i-1] * (2 * lambda * theta * 1 / delta) * Z1[i-1];
  }
  return p;
}
