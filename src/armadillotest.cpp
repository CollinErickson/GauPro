#include <RcppArmadillo.h>
using namespace arma;


// [[Rcpp::export]]
arma::mat cholC (arma::mat x) {
  return(chol(x)) ;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
*/
