#include <Rcpp.h>
using namespace Rcpp;

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
NumericMatrix Gaussian_hessianCC(NumericVector XX, NumericMatrix X, NumericVector Z, NumericMatrix Kinv, double mu_hat, NumericVector theta) {
  int n = X.nrow(); // # number of points already in design
  int d = XX.length(); // # input dimensions
  NumericVector Kinv_Zmu = Kinv * (Z - mu_hat); //#solve(R, Z - mu_hat)
  NumericMatrix d2ZZ(d, d);
  NumericVector exp_sum(n); // = numeric(n)
  for (int j = 0; j < n; j++ ){ //]j in 1:n) {
    exp_sum[j] = exp(-sum(theta * pow(XX - X(j, _), 2)));
  }
  for (int i=0; i<d; i++) { //(i in 1:d) { # diagonal points
    NumericVector d2K_dxidxi(n);// <- numeric(n)
    for (int j = 1; j < n; j++ ){ //(j in 1:n) {
      //d2Kj_dxidxi <- (-2 * theta[i] + 4 * theta[i]^2 * (XX[i] - X[j, i])^2) * exp_sum[j]
      //d2K_dxidxi[j] <- d2Kj_dxidxi
      d2K_dxidxi[j] = (-2 * theta[i] + 4 * pow(theta[i], 2) * pow(XX[i] - X[j, i], 2)) * exp_sum[j];
    }

    //tval <- t(d2K_dxidxi) %*% Kinv_Zmu
    //  d2ZZ[i, i] <- tval
    //d2ZZ[i, i] = d2K_dxidxi* Kinv_Zmu;
    for (int h=0; h<n; h++) {
      d2ZZ[i, i] += d2K_dxidxi(h) * Kinv_Zmu(h);
    }
    Rcout << i << d2ZZ[i, i] << std::endl;
  }
  if (d > 1) {
    for (int i=0; i < d-1; i++ ){ //(i in 1:(d-1)) { # off diagonal points
      for (int k=i+1; k < d; k++ ){ //(k in (i+1):d) {

        NumericVector d2K_dxidxk(n); // <- numeric(n)
        for (int j=0; j < n; j++ ){ //(j in 1:n) {

          //d2Kj_dxidxk = 4 * theta[i] * theta[k] * (XX[i] - X[j, i]) * (XX[k] - X[j, k]) * exp_sum[j];
          //d2K_dxidxk[j] = d2Kj_dxidxk;
          d2K_dxidxk[j] = 4 * theta[i] * theta[k] * (XX[i] - X[j, i]) * (XX[k] - X[j, k]) * exp_sum[j];
        }

        double tval = 0;
        for (int h=0; h<n; h++) {
          tval += d2K_dxidxk(h) * Kinv_Zmu(h);
        }
        //print("in C here");
        Rcout << i << k << tval << std::endl;
        //d2ZZ[i, k] = tval;
        //d2ZZ[k, i] = tval;
      }
    }
  }


  //hess <- d2ZZ
  return d2ZZ;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
#timesTwo(42)
*/
