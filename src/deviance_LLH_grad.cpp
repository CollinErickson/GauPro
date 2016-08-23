//#include <Rcpp.h>
#include <RcppArmadillo.h>
//  using namespace Rcpp;
using namespace arma;




// [[Rcpp::export]]
arma::vec deviance_LLH_grad_theta(arma::mat X, arma::mat Z, arma::mat K) {
  int N = X.n_rows;
  int D = X.n_cols;
  arma::vec grad = zeros(D);

  arma::mat Kchol = chol(K);
  double mu_hat_top = sum(sum(solve(trimatu(Kchol), solve(trimatl(Kchol.t()), Z))));
  arma::vec mu_hat_bottom_half = solve(trimatl(Kchol.t()), arma::ones(N));
  double mu_hat_bottom = sum(mu_hat_bottom_half.t() * mu_hat_bottom_half);
  double mu_hat = mu_hat_top / mu_hat_bottom;
  arma::vec y = Z - mu_hat;
  arma::vec Kinv_y = solve(trimatu(Kchol), solve(trimatl(Kchol.t()), Z - mu_hat));
  //double tmat = (sum(trans(Z - mu_hat) * (Kinv_y)));
  //double logdetK = 2 * sum(log(diagvec(Kchol)));
  //fngr[0] = logdetK + tmat;
  // Now calculate gradient
  double t1;
  arma::mat dK = zeros(N,N);
  // over thetas
  arma::mat Kcholtinv_dK = zeros(N, N);
  for(int i = 0; i < D; i++) {
    dK = K;
    for(int j = 0; j < N; j++) {
      for(int k = 0; k < N; k++) {
        // Doesn't like pow for some reason, but it works
        dK(j, k) = - pow(X(j,i) - X(k,i), 2) * dK(j, k);
      }
    }
    t1 = trace(Kinv_y * Kinv_y.t() * dK - solve(trimatu(Kchol), solve(trimatl(Kchol.t()), dK)));
    grad[i] = -2 * t1;
  }
  return grad;
}




// [[Rcpp::export]]
double deviance_LLH_grad_nug(arma::mat X, arma::mat Z, arma::mat K) {
  int N = X.n_rows;
  //int D = X.n_cols;
  //arma::vec fngr = zeros(1 + 1);

  arma::mat Kchol = chol(K);
  double mu_hat_top = sum(sum(solve(trimatu(Kchol), solve(trimatl(Kchol.t()), Z))));
  arma::vec mu_hat_bottom_half = solve(trimatl(Kchol.t()), arma::ones(N));
  double mu_hat_bottom = sum(mu_hat_bottom_half.t() * mu_hat_bottom_half);
  double mu_hat = mu_hat_top / mu_hat_bottom;
  arma::vec y = Z - mu_hat;
  arma::vec Kinv_y = solve(trimatu(Kchol), solve(trimatl(Kchol.t()), Z - mu_hat));
  //double tmat = (sum(trans(Z - mu_hat) * (Kinv_y)));
  //double logdetK = 2 * sum(log(diagvec(Kchol)));
  //fngr[0] = logdetK + tmat;
  // Now calculate gradient
  double t1;
  // for nugget
  t1 = 0; // just the trace
  arma::mat Kchol_inv = inv(trimatu(Kchol));
  for (int ii = 0; ii < N; ii++) {
    t1 += pow(Kinv_y(ii), 2) - dot(Kchol_inv.row(ii), Kchol_inv.row(ii));
  }
  return -t1;
  //fngr[1] = -t1;
  //return fngr;
}






// [[Rcpp::export]]
arma::vec deviance_LLH_grad_joint(arma::mat X, arma::mat Z, arma::mat K) {
  int N = X.n_rows;
  int D = X.n_cols;
  arma::vec grad = zeros(D + 1);

  arma::mat Kchol = chol(K);
  double mu_hat_top = sum(sum(solve(trimatu(Kchol), solve(trimatl(Kchol.t()), Z))));
  arma::vec mu_hat_bottom_half = solve(trimatl(Kchol.t()), arma::ones(N));
  double mu_hat_bottom = sum(mu_hat_bottom_half.t() * mu_hat_bottom_half);
  double mu_hat = mu_hat_top / mu_hat_bottom;
  arma::vec y = Z - mu_hat;
  arma::vec Kinv_y = solve(trimatu(Kchol), solve(trimatl(Kchol.t()), Z - mu_hat));
  //double tmat = (sum(trans(Z - mu_hat) * (Kinv_y)));
  //double logdetK = 2 * sum(log(diagvec(Kchol)));
  //fngr[0] = logdetK + tmat;
  // Now calculate gradient
  double t1;
  arma::mat dK = zeros(N,N);
  // over thetas
  arma::mat Kcholtinv_dK = zeros(N, N);
  for(int i = 0; i < D; i++) {
    dK = K;
    for(int j = 0; j < N; j++) {
      for(int k = 0; k < N; k++) {
        // Doesn't like pow for some reason, but it works
        dK(j, k) = - pow(X(j,i) - X(k,i), 2) * dK(j, k);
      }
    }
    t1 = trace(Kinv_y * Kinv_y.t() * dK - solve(trimatu(Kchol), solve(trimatl(Kchol.t()), dK)));
    grad[i] = -2 * t1;
  }
  // for nugget
  t1 = 0; // just the trace
  arma::mat Kchol_inv = inv(trimatu(Kchol));
  for (int ii = 0; ii < N; ii++) {
    t1 += pow(Kinv_y(ii), 2) - dot(Kchol_inv.row(ii), Kchol_inv.row(ii));
  }
  grad[D] = -t1;
  return grad;
}




// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
*/
