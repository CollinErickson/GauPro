//#include <Rcpp.h>
#include <RcppArmadillo.h>
//using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::vec deviance_grad_theta(arma::mat X, arma::mat K, arma::mat Kinv, arma::vec y) {
  int N = X.n_rows;
  int D = X.n_cols;
  arma::mat Kinv_y = Kinv * y;
  arma::mat t2amat = (trans(y) * Kinv_y);
  double t2a = -N / t2amat(0,0);
  arma::vec dD = zeros(D);
  double t1;
  arma::mat t2 = zeros(1,1);
  arma::mat dK = zeros(N,N);
  for(int i = 0; i < D; i++) {
    dK = K;
    for(int j = 0; j < N; j++) {
      for(int k = 0; k < N; k++) {
        // Doesn't like pow for some reason, but it works
        dK(j, k) = - pow(X(j,i) - X(k,i), 2) * dK(j, k);
      }
    }
    t1 = 0;
    for (int ii = 0; ii < N; ii++) {
      t1 += sum(Kinv.row(ii) * dK.col(ii));
    }
    t2 = t2a * trans(Kinv_y) * dK * Kinv_y;
    dD[i] = 2 * (t1 + t2(0,0));
  }
  return dD;
}

// [[Rcpp::export]]
double deviance_grad_nug(arma::mat X, arma::mat K, arma::mat Kinv, arma::vec y) {
  int N = X.n_rows;
  arma::mat Kinv_y = Kinv * y;
  arma::mat t2amat = (trans(y) * Kinv_y);
  double t2a = -N / t2amat(0,0);
  double dD = 0;
  double t1;
  arma::mat t2 = zeros(1,1);
  //arma::mat dK = diagmat(ones(N)); # dK is identity
  t1 = 0; // just the trace
  for (int ii = 0; ii < N; ii++) {
    t1 += Kinv(ii,ii); //sum(Kinv.row(ii) * dK.col(ii));
  }
  t2 = t2a * trans(Kinv_y) * Kinv_y;
  dD = t1 + t2(0,0);
  return dD;
}


// [[Rcpp::export]]
arma::vec deviance_grad_joint(arma::mat X, arma::mat K, arma::mat Kinv, arma::vec y) {
  int N = X.n_rows;
  int D = X.n_cols;
  arma::mat Kinv_y = Kinv * y;
  arma::mat t2amat = (trans(y) * Kinv_y);
  double t2a = -N / t2amat(0,0);
  arma::vec dD = zeros(D + 1);
  double t1;
  arma::mat t2 = zeros(1,1);
  arma::mat dK = zeros(N,N);
  // over thetas
  for(int i = 0; i < D; i++) {
    dK = K;
    for(int j = 0; j < N; j++) {
      for(int k = 0; k < N; k++) {
        // Doesn't like pow for some reason, but it works
        dK(j, k) = - pow(X(j,i) - X(k,i), 2) * dK(j, k);
      }
    }
    t1 = 0;
    for (int ii = 0; ii < N; ii++) {
      t1 += sum(Kinv.row(ii) * dK.col(ii));
    }
    t2 = t2a * trans(Kinv_y) * dK * Kinv_y;
    dD[i] = 2 * (t1 + t2(0,0));
  }
  // for nugget
  t1 = 0; // just the trace
  for (int ii = 0; ii < N; ii++) {
    t1 += Kinv(ii,ii); //sum(Kinv.row(ii) * dK.col(ii));
  }
  t2 = t2a * trans(Kinv_y) * Kinv_y;
  dD[D] = t1 + t2(0,0);
  return dD;
}




// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
#timesTwo(42)
*/
