#include <RcppArmadillo.h>
using namespace Rcpp;


//' Calculate gradfunc in optimization to speed up.
//' NEEDS TO APERM dC_dparams
//' Doesn't need to be exported, should only be useful in functions.
//' @param dC_dparams Derivative matrix for covariance function wrt kernel parameters
//' @param Cinv Inverse of covariance matrix
//' @param Cinv_yminusmu Vector that is the inverse of C times y minus the mean.
//' @return Vector, one value for each parameter
//' @examples
//' gradfuncarray(array(dim=c(2,4,4), data=rnorm(32)), matrix(rnorm(16),4,4), rnorm(4))
//' @export
// [[Rcpp::export]]
arma::vec gradfuncarray(arma::cube dC_dparams, arma::mat Cinv, arma::vec Cinv_yminusmu) {
  int d1 = dC_dparams.n_rows;
  int d2 = dC_dparams.n_cols;
  int d3 = dC_dparams.n_slices;
  arma::vec out(d1);
  double t1;
  double t2;
  for (int i = 0; i < d1; i++) {
    t1 = 0;
    t2 = 0;
    for (int j = 0; j < d2; j++) {
      for (int k = 0; k < d3; k++) {
        t1 += Cinv(j, k) * dC_dparams(i, j, k);
        t2 += Cinv_yminusmu(j) * dC_dparams(i, j, k) * Cinv_yminusmu(k);
      }
    }
    out(i) = t1 - t2;
  }
  return out;
}
