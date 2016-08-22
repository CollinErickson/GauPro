//#include <Rcpp.h>
#include <RcppArmadillo.h>
//  using namespace Rcpp;
using namespace arma;


arma::vec deviance_LLH_grad_nug(arma::mat X, arma::mat Z, arma::mat K) {
  int N = X.n_rows;
  arma::mat Kchol = chol(K);
  double mu_hat_top = sum(sum(solve(trimatu(Kchol), solve(trimatl(Kchol.t()), Z))));
  arma::vec mu_hat_bottom_half = solve(trimatl(Kchol.t()), arma::ones(N));
  double mu_hat_bottom = sum(mu_hat_bottom_half.t() * mu_hat_bottom_half);
  double mu_hat = mu_hat_top / mu_hat_bottom;
  arma::vec Kinv_y = solve(trimatu(Kchol), solve(trimatl(Kchol.t()), Z - mu_hat));
  double tmat = sum(trans(Z - mu_hat) * (Kinv_y));
  double logdetK = 2 * sum(log(diagvec(Kchol)));
  arma::vec dD = zeros(2);
  dD[0] = logdetK + tmat;//(0,0);

  dD[1] = -trace(Kinv_y * trans(Kinv_y));
  return dD;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
*/
