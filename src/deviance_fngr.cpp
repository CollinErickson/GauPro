//#include <Rcpp.h>
#include <RcppArmadillo.h>
//using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::vec deviance_fngr_theta(arma::mat X, arma::vec Z, arma::mat K) {
  int N = X.n_rows;
  int D = X.n_cols;
  arma::vec fngr = zeros(1 + D);

  arma::mat Kchol = chol(K);
  double mu_hat_top = sum(sum(solve(trimatu(Kchol), solve(trimatl(Kchol.t()), Z))));
  arma::vec mu_hat_bottom_half = solve(trimatl(Kchol.t()), arma::ones(N));
  double mu_hat_bottom = sum(mu_hat_bottom_half.t() * mu_hat_bottom_half);
  double mu_hat = mu_hat_top / mu_hat_bottom;
  arma::vec y = Z - mu_hat;
  arma::vec Kinv_y = solve(trimatu(Kchol), solve(trimatl(Kchol.t()), Z - mu_hat));
  double tmat = N * log(sum(trans(Z - mu_hat) * (Kinv_y)));
  double logdetK = 2 * sum(log(diagvec(Kchol)));
  //return logdetK + tmat;//(0,0);
  fngr[0] = logdetK + tmat;//(0,0);

  //arma::mat Kinv_y = Kinv * y;
  arma::mat t2amat = (trans(y) * Kinv_y);
  double t2a = -N / t2amat(0,0);
  //arma::vec dD = zeros(D);
  double t1;
  arma::mat t2 = zeros(1,1);
  arma::mat dK = zeros(N,N);
  arma::mat Kcholtinv_dK = zeros(N, N);
  for(int i = 0; i < D; i++) {
    dK = K;
    for(int j = 0; j < N; j++) {
      for(int k = 0; k < N; k++) {
        // Doesn't like pow for some reason, but it works
        dK(j, k) = - pow(X(j,i) - X(k,i), 2) * dK(j, k);
      }
    }

    // need to invert or something here, instead trying double solve and trace, no longer loop
    /*t1 = 0;
    Kcholtinv_dK = solve(trimatl(Kchol.t()), dK);
    for (int ii = 0; ii < N; ii++) {
      //t1 += sum(Kinv.row(ii) * dK.col(ii));
      t1 += sum(Kchol.row(ii) * Kcholtinv_dK.col(ii));
    }*/
    t1 = trace(solve(trimatu(Kchol), solve(trimatl(Kchol.t()), dK)));

    t2 = t2a * trans(Kinv_y) * dK * Kinv_y;
    //dD[i] = t1 + t2(0,0);
    fngr[1 + i] = 2 * (t1 + t2(0,0));
  }
  //return dD;
  return fngr;
}

// [[Rcpp::export]]
arma::vec deviance_fngr_nug(arma::mat X, arma::vec Z, arma::mat K) {
  int N = X.n_rows;


  arma::vec fngr = zeros(1 + 1);

  arma::mat Kchol = chol(K);
  double mu_hat_top = sum(sum(solve(trimatu(Kchol), solve(trimatl(Kchol.t()), Z))));
  arma::vec mu_hat_bottom_half = solve(trimatl(Kchol.t()), arma::ones(N));
  double mu_hat_bottom = sum(mu_hat_bottom_half.t() * mu_hat_bottom_half);
  double mu_hat = mu_hat_top / mu_hat_bottom;
  arma::vec y = Z - mu_hat;
  arma::vec Kinv_y = solve(trimatu(Kchol), solve(trimatl(Kchol.t()), Z - mu_hat));
  double tmat = N * log(sum(trans(Z - mu_hat) * (Kinv_y)));
  double logdetK = 2 * sum(log(diagvec(Kchol)));
  //return logdetK + tmat;//(0,0);
  fngr[0] = logdetK + tmat;//(0,0);


  //arma::mat Kinv_y = Kinv * y;
  arma::mat t2amat = (trans(y) * Kinv_y);
  double t2a = -N / t2amat(0,0);
  //double dD = 0;
  double t1;
  arma::mat t2 = zeros(1,1);
  //arma::mat dK = diagmat(ones(N)); # dK is identity
  t1 = 0; // just the trace
  arma::mat Kchol_inv = inv(trimatu(Kchol));
  for (int ii = 0; ii < N; ii++) {
    //t1 += Kinv(ii,ii); //sum(Kinv.row(ii) * dK.col(ii));
    t1 += dot(Kchol_inv.row(ii), Kchol_inv.row(ii));
  }
  t2 = t2a * trans(Kinv_y) * Kinv_y;
  //dD = t1 + t2(0,0);
  fngr[1] = t1 + t2(0,0);
  //return dD;
  return fngr;
}


// [[Rcpp::export]]
arma::vec deviance_fngr_joint(arma::mat X, arma::mat Z, arma::mat K) {
  int N = X.n_rows;
  int D = X.n_cols;
  arma::vec fngr = zeros(1 + D + 1);

  arma::mat Kchol = chol(K);
  double mu_hat_top = sum(sum(solve(trimatu(Kchol), solve(trimatl(Kchol.t()), Z))));
  arma::vec mu_hat_bottom_half = solve(trimatl(Kchol.t()), arma::ones(N));
  double mu_hat_bottom = sum(mu_hat_bottom_half.t() * mu_hat_bottom_half);
  double mu_hat = mu_hat_top / mu_hat_bottom;
  arma::vec y = Z - mu_hat;
  arma::vec Kinv_y = solve(trimatu(Kchol), solve(trimatl(Kchol.t()), Z - mu_hat));
  double tmat = N * log(sum(trans(Z - mu_hat) * (Kinv_y)));
  double logdetK = 2 * sum(log(diagvec(Kchol)));
  //return logdetK + tmat;//(0,0);
  fngr[0] = logdetK + tmat;//(0,0);


  //arma::mat Kinv_y = Kinv * y;
  arma::mat t2amat = (trans(y) * Kinv_y);
  double t2a = -N / t2amat(0,0);
  //arma::vec dD = zeros(D + 1);
  double t1;
  arma::mat t2 = zeros(1,1);
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
    /*t1 = 0;
    Kcholtinv_dK = solve(trimatl(Kchol.t()), dK);
    for (int ii = 0; ii < N; ii++) {
      //t1 += sum(Kinv.row(ii) * dK.col(ii));
      t1 += sum(Kchol.row(ii) * Kcholtinv_dK.col(ii)); // THIS DOESN'T WORK RIGHT
    }*/
    t1 = trace(solve(trimatu(Kchol), solve(trimatl(Kchol.t()), dK)));

    t2 = t2a * trans(Kinv_y) * dK * Kinv_y;
    //dD[i] = t1 + t2(0,0);
    fngr[1 + i] = 2 * (t1 + t2(0,0));
  }


  // for nugget
  t1 = 0; // just the trace
  arma::mat Kchol_inv = inv(trimatu(Kchol));
  for (int ii = 0; ii < N; ii++) {
    //t1 += Kinv(ii,ii); //sum(Kinv.row(ii) * dK.col(ii));
    t1 += dot(Kchol_inv.row(ii), Kchol_inv.row(ii));
  }
  t2 = t2a * trans(Kinv_y) * Kinv_y;
  //dD = t1 + t2(0,0);
  fngr[1 + D] = t1 + t2(0,0);

  //return dD;
  return fngr;
}




// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
#timesTwo(42)
*/
