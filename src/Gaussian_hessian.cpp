#include <RcppArmadillo.h>
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

//' Gaussian hessian in C
//'
//' @param XX point to find Hessian at
//' @param X matrix of data points
//' @param Z matrix of output
//' @param Kinv inverse of correlation matrix
//' @param mu_hat mean estimate
//' @param theta correlation parameters
//'
//' @return Hessian matrix
//' @export
// [[Rcpp::export]]
arma::mat Gaussian_hessianCC(arma::vec XX, arma::mat X, arma::vec Z, arma::mat Kinv, double mu_hat, arma::vec theta) {
  int n = X.n_rows; // # number of points already in design
  int d = XX.n_elem; // # input dimensions
  //Rcpp::Rcout << "at 1" << std::endl;
  arma::vec Kinv_Zmu = Kinv * (Z - mu_hat); //#solve(R, Z - mu_hat)

  // Initialize with zeros to avoid errors, d2ZZ was actually giving error on valgrind
  //arma::mat d2ZZ(d, d);
  //arma::vec exp_sum(n); // = numeric(n)
  arma::mat d2ZZ = zeros(d, d);
  arma::vec exp_sum = zeros(n); // = numeric(n)

  for (int j = 0; j < n; j++ ){ //)j in 1:n) {
    //Rcpp::Rcout << "at 2" << std::endl;
    exp_sum(j) = exp(-sum(theta % pow(XX - X.row(j).t(), 2)));
  }
  //Rcpp::Rcout << "at 3" << std::endl;
  //arma::vec d2K_dxidxi(n);// <- numeric(n)
  arma::vec d2K_dxidxi = zeros(n);// Initialized with zeros to be safe
  for (int i=0; i<d; i++) { //(i in 1:d) { # diagonal points
    // d2K_dxidxi used to initialize here. It overwrites all so I don't need to set each element back to 0.
    for (int j = 0; j < n; j++ ){ //(j in 1:n) {
      //d2Kj_dxidxi <- (-2 * theta(i) + 4 * theta(i)^2 * (XX(i) - X(j, i))^2) * exp_sum(j)
      //d2K_dxidxi(j) <- d2Kj_dxidxi
      d2K_dxidxi(j) = (-2 * theta(i) + 4 * pow(theta(i), 2) * pow(XX(i) - X(j, i), 2)) * exp_sum(j);
    }

    //tval <- t(d2K_dxidxi) %*% Kinv_Zmu
    //  d2ZZ(i, i) <- tval
    //d2ZZ(i, i) = d2K_dxidxi* Kinv_Zmu;
    for (int h=0; h<n; h++) {
      d2ZZ(i, i) += d2K_dxidxi(h) * Kinv_Zmu(h);
    }
    //Rcout << i << d2ZZ(i, i) << std::endl;
  }
  //Rcout << i << Kinv_Zmu << std::endl;
  //Rcout << d2ZZ << std::endl;
  //arma::vec d2K_dxidxk(n); // initialize outside of if and for to avoid problems
  arma::vec d2K_dxidxk = zeros(n); // initialize with zeros
  double tval = 0; // Initialize out here so it doesn't reinitialize each iteration
  if (d > 1) {
    for (int i=0; i < d-1; i++ ){ //(i in 1:(d-1)) { # off diagonal points
      for (int k=i+1; k < d; k++ ){ //(k in (i+1):d) {

        //arma::vec d2K_dxidxk(n); // <- numeric(n)
        for (int j=0; j < n; j++ ){ //(j in 1:n) {

          //d2Kj_dxidxk = 4 * theta(i) * theta(k) * (XX(i) - X(j, i)) * (XX(k) - X(j, k)) * exp_sum(j);
          //d2K_dxidxk(j) = d2Kj_dxidxk;
          d2K_dxidxk(j) = 4 * theta(i) * theta(k) * (XX(i) - X(j, i)) * (XX(k) - X(j, k)) * exp_sum(j);
        }

        tval = 0;
        for (int h=0; h<n; h++) {
          tval += d2K_dxidxk(h) * Kinv_Zmu(h);
        }
        //print("in C here");
        //Rcout << i << k << tval << std::endl;
        d2ZZ(i, k) = tval;
        d2ZZ(k, i) = tval;
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
