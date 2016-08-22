//#include <Rcpp.h>
#include <RcppArmadillo.h>
//  using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
double deviance_LLH(arma::vec theta, double nug, arma::mat X, arma::mat Z, arma::mat K) {
  // Twice as fast to this compared to devianceC or just R version
  int N = X.n_rows;
  arma::mat Kchol = chol(K);
  double mu_hat_top = sum(sum(solve(trimatu(Kchol), solve(trimatl(Kchol.t()), Z))));
  arma::vec mu_hat_bottom_half = solve(trimatl(Kchol.t()), arma::ones(N));
  double mu_hat_bottom = sum(mu_hat_bottom_half.t() * mu_hat_bottom_half);
  double mu_hat = mu_hat_top / mu_hat_bottom;
  arma::vec Kinv_y = solve(trimatu(Kchol), solve(trimatl(Kchol.t()), Z - mu_hat));
  double tmat = sum(trans(Z - mu_hat) * (Kinv_y));
  double logdetK = 2 * sum(log(diagvec(Kchol)));
  return logdetK + tmat;//(0,0);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
#timesTwo(42)
*/
