//#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
//using namespace arma;


//' Correlation ordered factor  matrix in C (symmetric)
//' @param x Matrix x
//' @param theta Theta vector
//' @param xindex Index to use
//' @param offdiagequal What to set off-diagonal values with matching values to.
//' @return Correlation matrix
//' @export
//' @examples
//' corr_orderedfactor_matrix_symC(matrix(c(1,.5, 2,1.6, 1,0),ncol=2,byrow=TRUE),
//'                               c(1.5,1.8), 1, 1-1e-6)
//' corr_orderedfactor_matrix_symC(matrix(c(0,0,0,1,0,0,0,2,0,0,0,3,0,0,0,4),
//'                                      ncol=4, byrow=TRUE),
//'   c(0.101, -0.714, 0.114, -0.755, 0.117, -0.76, 0.116, -0.752),
//'   4, 1-1e-6) * 6.85
// [[Rcpp::export]]
NumericMatrix corr_orderedfactor_matrix_symC(NumericMatrix x,
                                             NumericVector theta,
                                             int xindex,
                                             double offdiagequal) {
  int nrow = x.nrow();
  //int nsum = x.ncol();
  NumericMatrix out(nrow, nrow);
  // int xindoffset, yindoffset;
  int xlev;
  int ylev;
  int a, b;
  double total;
  for (int i = 0; i < nrow - 1; i++) {
    for (int j = i + 1; j < nrow; j++) {
      xlev = x(i, xindex - 1);
      ylev = x(j, xindex - 1);
      if (xlev == ylev) {
        total = offdiagequal;
      } else {
        // xindoffset = (xlev - 1) * latentdim;
        // yindoffset = (ylev - 1) * latentdim;
        // int a = min(xlev, ylev);
        // int b = max(xlev, ylev) - 1;
        if (xlev <= ylev) {
          a = xlev - 1;
          b = ylev - 2;
        } else {
          a = ylev - 1;
          b = xlev - 2;
        }
        total = 0;
        for(int k = a; k <= b; ++k) {
          total += theta[k];
        }
        total = exp(-pow(total, 2));
      }
      out(i, j) = total;
      out(j, i) = total; // since symmetric
    }
  }
  // Diagonal
  for (int i = 0; i < nrow; i++) {
    out(i, i) = 1;
  }
  return out;
}

/*** R
6.85 * corr_orderedfactor_matrix_symC(
  matrix(c(0,0,0,1,0,0,0,2,0,0,0,3,0,0,0,4), ncol=4, byrow=T),
  c(1,1,2),
  4, 1-1e-6)
*/


//' Correlation ordered factor matrix in C (symmetric)
//' @param x Matrix x
//' @param y Matrix y
//' @param theta Theta vector
//' @param xindex Index to use
//' @param offdiagequal What to set off-diagonal values with matching values to.
//' @return Correlation matrix
//' @export
//' @examples
//' corr_orderedfactor_matrixmatrixC(matrix(c(1,.5, 2,1.6, 1,0),ncol=2,byrow=TRUE),
//'                                 matrix(c(2,1.6, 1,0),ncol=2,byrow=TRUE),
//'                                 c(1.5,1.8), 1, 1-1e-6)
//' corr_orderedfactor_matrixmatrixC(matrix(c(0,0,0,1,0,0,0,2,0,0,0,3,0,0,0,4),
//'                                   ncol=4, byrow=TRUE),
//'                                 matrix(c(0,0,0,2,0,0,0,4,0,0,0,1),
//'                                   ncol=4, byrow=TRUE),
//'   c(0.101, -0.714, 0.114, -0.755, 0.117, -0.76, 0.116, -0.752),
//'   4, 1-1e-6) * 6.85
// [[Rcpp::export]]
NumericMatrix corr_orderedfactor_matrixmatrixC(NumericMatrix x, NumericMatrix y,
                                               NumericVector theta,
                                               int xindex,
                                               double offdiagequal) {
  int nrow = x.nrow();
  int ncol = y.nrow();
  //int nsum = x.ncol();
  NumericMatrix out(nrow, ncol);
  // int xindoffset, yindoffset;
  int xlev;
  int ylev;
  int a, b;
  double total;
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      xlev = x(i, xindex - 1);
      ylev = y(j, xindex - 1);
      if (xlev == ylev) {
        total = offdiagequal;
      } else {
        // xindoffset = (xlev - 1) * latentdim;
        // yindoffset = (ylev - 1) * latentdim;
        if (xlev <= ylev) {
          a = xlev - 1;
          b = ylev - 2;
        } else {
          a = ylev - 1;
          b = xlev - 2;
        }
        total = 0;
        // double latx, laty;
        for(int k = a; k <= b; ++k) {
          // latx = theta[xindoffset + k];
          // laty = theta[yindoffset + k];
          // total += pow(latx - laty, 2);
          total += theta[k];
        }
        total = exp(-pow(total, 2));
      }
      out(i, j) = total;
    }
  }
  return out;
}


// Trying to get C_dC for Gaussian kernel
//' Derivative of covariance matrix of X with respect to kernel
//' parameters for the Ordered Factor Kernel
//' @param x Matrix x
//' @param pf pf vector
//' @param C_nonug cov mat without nugget
//' @param s2_est whether s2 is being estimated
//' @param p_est Whether theta/beta is being estimated
//' @param lenparams_D Number of parameters the derivative is being calculated for
//' @param s2_nug s2 times the nug
//' @param xindex Which column of x is the indexing variable
//' @param nlevels Number of levels
//' @param s2 Value of s2
//' @return Correlation matrix
//' @export
// [[Rcpp::export]]
arma::cube kernel_orderedFactor_dC(
    arma::mat x, arma::vec pf, arma::mat C_nonug,
    bool s2_est, bool p_est, int lenparams_D, double s2_nug,
    int xindex, int nlevels, double s2) {
  int nrow = x.n_rows;
  //int nsum = x.n_cols;
  arma::cube dC_dparams(lenparams_D, nrow, nrow);

  if (s2_est) {
    for (int i = 0; i < nrow - 1; i++) {
      for (int j = i + 1; j < nrow; j++) {
        dC_dparams(lenparams_D - 1,i,j) = C_nonug(i,j) * log(10.0);
        dC_dparams(lenparams_D - 1,j,i) = dC_dparams(lenparams_D - 1,i,j);
      }
      dC_dparams(lenparams_D - 1, i, i) = (C_nonug(i,i) + s2_nug) * log(10.0);
    }
    dC_dparams(lenparams_D - 1, nrow - 1, nrow - 1) = (
      C_nonug(nrow - 1, nrow - 1) + s2_nug) * log(10.0);
  }
  if (p_est) {
    int xlev, ylev;
    int a, b;
    for (int i = 0; i < nrow - 1; i++) {
      for (int j = i + 1; j < nrow; j++) {
        xlev = x(i, xindex);
        ylev = x(j, xindex);
        if (xlev == ylev) {
          // Params don't affect
        } else {
          if (xlev <= ylev) {
            a = xlev - 1;
            b = ylev - 2;
          } else {
            a = ylev - 1;
            b = xlev - 2;
          }
          for (int k = 0; k < nlevels - 1; k++) {
            if (a <= k && k <= b) {
              // dC_dparams(k_ind,i,j) = -2 * out * (latentx_l - latenty_l);
              dC_dparams(k,i,j) = (
                -2 * sqrt(-log(C_nonug(i, j) / s2)) * C_nonug(i, j));
              dC_dparams(k,j,i) = dC_dparams(k,i,j);
            }
            // if (xlev>1.5 && xlev==k && ylev !=k) {
            //   double p_dist2 = 0;
            //   for (int l=0; l<latentdim; l++) {
            //     double latentx_l = pf[(xlev-1)*latentdim + l];
            //     double latenty_l = pf[(ylev-1)*latentdim + l];
            //     p_dist2 += std::pow(latentx_l - latenty_l, 2);
            //   }
            //   // s2 or s2_nug here?
            //   double out = s2 * exp(-p_dist2);
            //
            //   for (int l=0; l<latentdim; l++) {
            //     int k_ind = (k-2)*latentdim + l;
            //     double latentx_l = pf[(xlev-1)*latentdim + l];
            //     double latenty_l = pf[(ylev-1)*latentdim + l];
            //     dC_dparams(k_ind,i,j) = -2 * out * (latentx_l - latenty_l);
            //     dC_dparams(k_ind,j,i) = dC_dparams(k_ind,i,j);
            //   }
            // }
          }
        }
      }
    }
    for (int k = 0; k < lenparams_D - 1; k++) {
      for (int i = 0; i < nrow; i++) { //# Get diagonal set to zero
        dC_dparams(k,i,i) = 0;
      }
    }
  }

  return dC_dparams;
}
