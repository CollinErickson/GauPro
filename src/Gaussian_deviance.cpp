//  #include <Rcpp.h>
#include <RcppArmadillo.h>
#include "corr.h"
//using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
double Gaussian_deviance_part(arma::vec theta, double nug, arma::mat X, arma::mat Z, arma::mat Kinv) {
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


//' Calculate the Gaussian deviance in C
//' @param X Matrix X
//' @param Z Matrix Z
//' @param theta Theta vector
//' @param nug Nugget
//' @return Correlation matrix
//' @examples
//' Gaussian_devianceC(c(1,1), 1e-8, matrix(c(1,0,0,1),2,2), matrix(c(1,0),2,1))
//' @export
// [[Rcpp::export]]
double Gaussian_devianceC(arma::vec theta, double nug, arma::mat X, arma::mat Z) {
  // Twice as fast to this compared to devianceC or just R version
  int N = X.n_rows;

  arma::vec nug_vec(N);
  for (int i =0; i<N; i++) {nug_vec[i] = nug;};
  arma::mat K = (corr_gauss_matrix_sym_armaC(X, theta)) + arma::diagmat(nug_vec);//  diag(nug, self$N);

  arma::mat Kchol = chol(K);
  double mu_hat_top = sum(sum(solve(trimatu(Kchol), solve(trimatl(Kchol.t()), Z))));
  arma::vec mu_hat_bottom_half = solve(trimatl(Kchol.t()), arma::ones(N));
  double mu_hat_bottom = sum(mu_hat_bottom_half.t() * mu_hat_bottom_half);
  double mu_hat = mu_hat_top / mu_hat_bottom;
  arma::vec Kinv_y = solve(trimatu(Kchol), solve(trimatl(Kchol.t()), Z - mu_hat));
  double tmat = N * log(sum(trans(Z - mu_hat) * (Kinv_y)));
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
