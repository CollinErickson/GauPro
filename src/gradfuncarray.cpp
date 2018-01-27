//This ended up being 5x slower than using apply and R matrix operations,
//  so not going to use it
// Below were in GauPro_init.c

//extern SEXP _GauPro_gradfuncarray(SEXP, SEXP, SEXP);
//extern SEXP _GauPro_gradfuncarray2(SEXP, SEXP, SEXP);

//{"_GauPro_gradfuncarray",               (DL_FUNC) &_GauPro_gradfuncarray,         3},
//{"_GauPro_gradfuncarray2",               (DL_FUNC) &_GauPro_gradfuncarray2,         3},


#include <RcppArmadillo.h>
using namespace Rcpp;


//' Calculate gradfunc in optimization to speed up.
//' NEEDS TO APERM dC_dparams
//' Doesn't need to be exported, should only be useful in functions.
//' @param dC_dparams Derivative matrix for covariance function wrt kernel parameters
//' @param C Covariance matrix
//' @param Cinv_yminusmu Vector that is the inverse of C times y minus the mean.
//' @return Vector, one value for each parameter
//' @examples
//' # corr_gauss_dCdX(matrix(c(1,0,0,1),2,2),c(1,1))
//' @export
// [[Rcpp::export]]
arma::vec gradfuncarray(arma::cube dC_dparams, arma::mat Cinv, arma::vec Cinv_yminusmu) {
  int d1 = dC_dparams.n_rows;
  int d2 = dC_dparams.n_cols;
  int d3 = dC_dparams.n_slices;
  arma::vec out(d1);
  double t1;
  double t2;
  // arma::mat tmat;
  // arma::mat tmat2;
  // arma::vec tvec1;
  for (int i = 0; i < d1; i++) {

    // t1 <- sum(Cinv * t(di))
    // t2 <- sum(Cinv_yminusmu * (di %*% Cinv_yminusmu))
    // t1 - t2
    t1 = 0;
    t2 = 0;
    for (int j = 0; j < d2; j++) {
      for (int k = 0; k < d3; k++) {
        t1 += Cinv(j, k) * dC_dparams(i, j, k);
        t2 += Cinv_yminusmu(j) * dC_dparams(i, j, k) * Cinv_yminusmu(k);
      }
    }

    // tmat = dC_dparams.subcube(i,0,0,i,d2-1,d3-1);
    // t1 = arma::sum((arma::solve(C, dC_dparams.subcube(0,0,i,d1-1,d2-2,i))).diag());
    // t2 = arma::sum(Cinv_yminusmu, (dC_dparams.subcube(0,0,i,d1-1,d2-2,i) * Cinv_yminusmu));
    // tmat2 = arma::solve(C, dC_dparams.slice(i));
    // t1 = arma::sum(tmat2.diag());
    // tvec1 = dC_dparams.slice(i) * Cinv_yminusmu;
    // t2 = arma::dot(Cinv_yminusmu, tvec1); // * (tmat * Cinv_yminusmu));
    out(i) = t1 - t2;
  }

  return out;
}

