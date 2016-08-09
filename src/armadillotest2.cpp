#include <RcppArmadillo.h>
using namespace arma;


// [[Rcpp::export]]
arma::mat solveC (arma::mat A, arma::vec b) {
  return(solve(A,b)) ;
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
*/
