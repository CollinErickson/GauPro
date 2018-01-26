#include <RcppArmadillo.h>
using namespace Rcpp;

//' Correlation Gaussian matrix gradient in C using Armadillo
//' @param XX Matrix XX to get gradient for
//' @param X Matrix X GP was fit to
//' @param theta Theta vector
//' @param s2 Variance parameter
//' @return 3-dim array of correlation derivative
//' @examples
//' # corr_gauss_dCdX(matrix(c(1,0,0,1),2,2),c(1,1))
//' @export
// [[Rcpp::export]]
arma::cube corr_gauss_dCdX(arma::mat XX, arma::mat X, arma::vec theta, double s2) {
  int nn = XX.n_rows;
  int d = XX.n_cols;
  int n = X.n_rows;
  arma::cube dC_dx(nn, d, n);
  double tsum = 0;
  for (int i = 0; i < nn; i++) {
    for (int j = 0; j < d; j++) {
      for (int k = 0; k < n; k++) {
        // dC_dx(i, j, k) = -2 * theta(j) * (XX(i, j) - X(k, j)) * s2 * exp(-sum(theta * (XX(i,) - X(k,)) ^ 2));
        tsum = 0;
        for (int l = 0; l < d; l++) {
          // sum(theta * (XX(i,) - X(k,)) ^ 2)
          tsum += theta(l) * pow(XX(i,l) - X(k,l), 2);
        }
        dC_dx(i, j, k) = -2 * theta(j) * (XX(i, j) - X(k, j)) * s2 * exp(-tsum);
      }
    }
  }

  return dC_dx;
}

