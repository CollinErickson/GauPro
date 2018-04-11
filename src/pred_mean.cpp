
#include <RcppArmadillo.h>
//  using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
arma::vec pred_meanC(arma::mat XX, arma::mat kx_xx, double mu_hat, arma::mat Kinv, arma::mat Z) {
  return mu_hat + trans(kx_xx) * Kinv * (Z - mu_hat);
}

// [[Rcpp::export]]
arma::vec pred_var(arma::mat XX, arma::mat kxx, arma::mat kx_xx, double s2_hat, arma::mat Kinv, arma::mat Z) {
  return s2_hat * arma::diagvec(kxx - trans(kx_xx) * Kinv * kx_xx);
}

// [[Rcpp::export]]
arma::mat pred_cov(arma::mat XX, arma::mat kxx, arma::mat kx_xx, double s2_hat, arma::mat Kinv, arma::mat Z) {
  return s2_hat * (kxx - trans(kx_xx) * Kinv * kx_xx);
}


// Creating versions that take in mu_hat as vector

// Does prediction using full Kinv, takes O(n^2)
// [[Rcpp::export]]
arma::vec pred_meanC_mumat(arma::mat XX, arma::mat kx_xx, arma::mat mu_hatX, arma::mat mu_hatXX, arma::mat Kinv, arma::mat Z) {
  return mu_hatXX + trans(kx_xx) * Kinv * (Z - mu_hatX);
}

// Now I precalculate Kinv_Z_minus_mu_hatX since it doesn't depend on XX.
// Makes it O(n) instead of O(n^2)
// [[Rcpp::export]]
arma::vec pred_meanC_mumat_fast(arma::mat XX, arma::mat kx_xx, arma::vec Kinv_Z_minus_mu_hatX, arma::mat mu_hatXX) {
  // return mu_hatXX + sum(kx_xx * Kinv_Z_minus_mu_hatX);
  arma::vec tvec = zeros(mu_hatXX.n_elem);
  int ncols = kx_xx.n_cols;
  int nrows = kx_xx.n_rows;
  for (int i = 0; i < ncols; i++) {
    tvec(i) = mu_hatXX(i, 0);
    for (int j = 0; j < nrows; j++) {
      tvec(i) += kx_xx(j, i) * Kinv_Z_minus_mu_hatX(j);
    }
  }
  return tvec;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
#timesTwo(42)
*/
