// #include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;


//' Cube multiply over first dimension
//'
//' The result is transposed since that is what apply will give you
//'
//' @param cub A cube (3D array)
//' @param v A vector
//' @return Transpose of multiplication over first dimension of cub time v
//' @examples
//' d1 <- 10
//' d2 <- 1e2
//' d3 <- 2e2
//' aa <- array(data = rnorm(d1*d2*d3), dim = c(d1, d2, d3))
//' bb <- rnorm(d3)
//' t1 <- apply(aa, 1, function(U) {U%*%bb})
//' t2 <- arma_mult_cube_vec(aa, bb)
//' dd <- t1 - t2
//'
//' summary(dd)
//' image(dd)
//' table(dd)
//' # microbenchmark::microbenchmark(apply(aa, 1, function(U) {U%*%bb}),
//' #                                arma_mult_cube_vec(aa, bb))
//' @export
// [[Rcpp::export]]
arma::mat arma_mult_cube_vec(arma::cube cub, arma::vec v) {
  int d1 = cub.n_rows;
  int d2 = cub.n_cols;
  int d3 = cub.n_slices; // equals v.n_rows
  // int d4 = v.n_cols // should be 1
  arma::mat out(d2, d1); // transposed version
  double total;
  for (int i = 0; i < d1; i++) {
    for (int j = 0; j < d2; j++) {
      total = 0;
      for(int k = 0; k < d3; ++k) {
        total += cub(i,j,k) * v(k);
      }

      // out(i, j) = total;
      out(j, i) = total; // Transposed version
    }
  }
  return out;
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
*/
