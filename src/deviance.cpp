//#include <Rcpp.h>
#include <RcppArmadillo.h>
//  using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
double deviance_part(arma::vec theta, double nug, arma::mat X, arma::mat Z, arma::mat Kinv) {
  // Not faster than using R, no need for this
  int N = X.n_rows;
  //double sumKinv = sum(sum(Kinv));
  double mu_hat = sum(sum(Kinv * Z)) / sum(sum(Kinv)) ;
  //arma::mat s2_hat_mat = trans(Z - mu_hat) * Kinv * (Z - mu_hat) / N;
  //double s2_hat = s2_hat_mat(1,1);
  //double logdetK = 2 * sum(log(diag(Kchol)));
  arma::mat tmat = N * log(trans(Z - mu_hat) * (Kinv * (Z - mu_hat)));
  return tmat(0,0);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
#timesTwo(42)
*/
