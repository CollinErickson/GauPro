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


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
#timesTwo(42)
*/
