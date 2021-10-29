# CE 10/26/21
# gradfuncarray using Rcpp kept crashing RStudio on my laptop.
# Just going to use R code and hope it doesn't crash.
# This probably isn't a timing bottleneck anyways.

#' Calculate gradfunc in optimization to speed up.
#' NEEDS TO APERM dC_dparams
#' Doesn't need to be exported, should only be useful in functions.
#' @param dC_dparams Derivative matrix for covariance function wrt kernel parameters
#' @param Cinv Inverse of covariance matrix
#' @param Cinv_yminusmu Vector that is the inverse of C times y minus the mean.
#' @return Vector, one value for each parameter
#' @examples
#' a1 <- array(dim=c(2,4,4), data=rnorm(32))
#' a2 <- matrix(rnorm(16),4,4)
#' a3 <- rnorm(4)
#' #gradfuncarray(a1, a2, a3)
#' #gradfuncarrayR(a1, a2, a3)
#' @export
gradfuncarrayR <- function(dC_dparams, Cinv, Cinv_yminusmu) {
  # Rcout << dC_dparams;
  # Rcout << "\n\nCinv\n";
  # Rcout << Cinv;
  # Rcout << "\n\nCinv_yminusmu\n";
  # Rcout << Cinv_yminusmu;
  # int d1 = dC_dparams.n_rows;
  # int d2 = dC_dparams.n_cols;
  # int d3 = dC_dparams.n_slices;
  # arma::vec out(d1);
  # double t1;
  # double t2;
  d1 = dim(dC_dparams)[1]
  d2 = dim(dC_dparams)[2]
  d3 = dim(dC_dparams)[3]
  out <- numeric(d1)
  if (d1 == 0L) {
    return(out)
  }
  for (i in 1:d1) { #int i = 0; i < d1; i++) {
    t1 = 0;
    t2 = 0;
    for (j in 1:d2) { #int j = 0; j < d2; j++) {
      for (k in 1:d3) { #int k = 0; k < d3; k++) {
        # t1 += Cinv(j, k) * dC_dparams(i, j, k);
        # t2 += Cinv_yminusmu(j) * dC_dparams(i, j, k) * Cinv_yminusmu(k);
        t1 <- t1 + Cinv[j, k] * dC_dparams[i, j, k];
        t2 <- t2 + Cinv_yminusmu[j] * dC_dparams[i, j, k] * Cinv_yminusmu[k];
      }
    }
    out[i] = t1 - t2;
  }
  return( out);
}

