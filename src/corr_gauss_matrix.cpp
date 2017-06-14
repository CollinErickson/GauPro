//#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
//using namespace arma;

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
NumericMatrix corr_gauss_matrixC(NumericMatrix x, NumericMatrix y, NumericVector theta) {
  int nrow = x.nrow(), ncol = y.nrow();
  int nsum = x.ncol();
  NumericMatrix out(nrow, ncol);

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {

      double total = 0;
      for(int k = 0; k < nsum; ++k) {
        total += theta[k] * pow((x(i,k) - y(j,k)), 2.0);
      }
      total = exp(-total);

      out(i, j) = total;
    }
  }
  return out;
}

//' Correlation Gaussian matrix in C (symmetric)
//' @param x Matrix x
//' @param theta Theta vector
//' @return Correlation matrix
//' @export
//' @examples
//' corr_gauss_matrix_symC(matrix(c(1,0,0,1),2,2),c(1,1))
// [[Rcpp::export]]
NumericMatrix corr_gauss_matrix_symC(NumericMatrix x, NumericVector theta) {
  int nrow = x.nrow();
  int nsum = x.ncol();
  NumericMatrix out(nrow, nrow);

  for (int i = 0; i < nrow - 1; i++) {
    for (int j = i + 1; j < nrow; j++) {

      double total = 0;
      for(int k = 0; k < nsum; ++k) {
        total += theta[k] * pow((x(i,k) - x(j,k)), 2.0);
      }
      total = exp(-total);

      out(i, j) = total;
      out(j, i) = total; // since symmetric
    }
  }
  for (int i = 0; i < nrow; i++) {
    out(i, i) = 1;
  }
  return out;
}




// [[Rcpp::export]]
NumericVector corr_gauss_matrixvecC(NumericMatrix x, NumericVector y, NumericVector theta) {
  int nrow = x.nrow(); //, ncol = y.nrow();
  int nsum = x.ncol();
  NumericVector out(nrow);

  for (int i = 0; i < nrow; i++) {
    double total = 0;
    for(int k = 0; k < nsum; ++k) {
      total += theta[k] * pow((x(i,k) - y(k)), 2.0);
    }
    total = exp(-total);

    out(i) = total;
  }
  return out;
}



//' Correlation Gaussian matrix in C using Armadillo (symmetric)
//' @param x Matrix x
//' @param theta Theta vector
//' @return Correlation matrix
//' @examples
//' corr_gauss_matrix_sym_armaC(matrix(c(1,0,0,1),2,2),c(1,1))
//' @export
// [[Rcpp::export]]
arma::mat corr_gauss_matrix_sym_armaC(arma::mat x, arma::vec theta) {
  int nrow = x.n_rows;
  int nsum = x.n_cols;
  arma::mat out(nrow, nrow);

  for (int i = 0; i < nrow - 1; i++) {
    for (int j = i + 1; j < nrow; j++) {

      double total = 0;
      for(int k = 0; k < nsum; ++k) {
        total += theta[k] * pow((x(i,k) - x(j,k)), 2.0);
      }
      total = exp(-total);

      out(i, j) = total;
      out(j, i) = total;
    }
  }
  for (int i = 0; i < nrow; i++) {
    out(i, i) = 1;
  }
  return out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
#timesTwo(42)
*/
