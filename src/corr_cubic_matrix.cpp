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
NumericMatrix corr_cubic_matrixC(NumericMatrix x, NumericMatrix y, NumericVector theta) {
  int nrow = x.nrow(), ncol = y.nrow();
  int nsum = x.ncol();
  NumericMatrix out(nrow, ncol);

  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {

      double total = 1;
      for(int k = 0; k < nsum; ++k) {
        // total += theta[k] * pow((x(i,k) - y(j,k)), 2.0);
        double d = fabs(x(i,k) - y(j,k)) / theta[k];
        double r = 0;
        if (d <= .5) {
          r = 1-6*pow(d, 2.0)+6*pow(d, 3.0);
        } else if (d <= 1) {
          r = 2*pow(1-d, 3.0);
        } else {
          r = 0;
        }
        total *= r;
      }
      out(i, j) = total;
    }
  }
  return out;
}

//' Correlation Cubic matrix in C (symmetric)
//' @param x Matrix x
//' @param theta Theta vector
//' @return Correlation matrix
//' @export
//' @examples
//' corr_cubic_matrix_symC(matrix(c(1,0,0,1),2,2),c(1,1))
// [[Rcpp::export]]
NumericMatrix corr_cubic_matrix_symC(NumericMatrix x, NumericVector theta) {
   int nrow = x.nrow();
   int nsum = x.ncol();
   NumericMatrix out(nrow, nrow);

   for (int i = 0; i < nrow - 1; i++) {
     for (int j = i + 1; j < nrow; j++) {

       double total = 1;
       for(int k = 0; k < nsum; ++k) {
         // total += theta[k] * pow((x(i,k) - x(j,k)), 2.0);
         double d = fabs(x(i,k) - x(j,k)) / theta[k];
         double r = 0;
         if (d <= .5) {
           r = 1-6*pow(d, 2.0) + 6*pow(d, 3.0);
         } else if (d <= 1) {
           r = 2*pow(1-d, 3.0);
         } else {
           r = 0;
         }
         total *= r;
       }
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
NumericVector corr_cubic_matrixvecC(NumericMatrix x, NumericVector y,
                                    NumericVector theta) {
  int nrow = x.nrow(); //, ncol = y.nrow();
  int nsum = x.ncol();
  NumericVector out(nrow);

  for (int i = 0; i < nrow; i++) {
    double total = 1;
    for(int k = 0; k < nsum; ++k) {
      // total += theta[k] * pow((x(i,k) - y(k)), 2.0);
      double d = fabs(x(i,k) - y(k)) / theta[k];
      double r = 0;
      if (d <= .5) {
        r = 1-6*pow(d, 2.0)+6*pow(d, 3.0);
      } else if (d <= 1) {
        r = 2*pow(1-d, 3.0);
      } else {
        r = 0;
      }
      total *= r;
    }
    out(i) = total;
  }
  return out;
}



//' Derivative of cubic kernel covariance matrix in C
//' @param x Matrix x
//' @param theta Theta vector
//' @param C_nonug cov mat without nugget
//' @param s2_est whether s2 is being estimated
//' @param beta_est Whether theta/beta is being estimated
//' @param lenparams_D Number of parameters the derivative is being calculated for
//' @param s2_nug s2 times the nug
//' @param s2 s2
//' @return Correlation matrix
//' @export
// [[Rcpp::export]]
arma::cube kernel_cubic_dC(arma::mat x, arma::vec theta, arma::mat C_nonug,
                            bool s2_est, bool beta_est, int lenparams_D,
                            double s2_nug, double s2) {
   int nrow = x.n_rows;
   int nsum = x.n_cols;
   arma::cube dC_dparams(lenparams_D, nrow, nrow);
   double log10 = log(10.0);

   if (s2_est) {
     // dC_dparams(lenparams_D,,) = C * log10;
     for (int i = 0; i < nrow - 1; i++) {
       for (int j = i + 1; j < nrow; j++) {
         dC_dparams(lenparams_D - 1,i,j) = C_nonug(i,j) * log10;
         dC_dparams(lenparams_D - 1,j,i) = dC_dparams(lenparams_D - 1,i,j);
       }
       dC_dparams(lenparams_D - 1, i, i) = (C_nonug(i,i) + s2_nug) * log10;
     }
     dC_dparams(lenparams_D - 1, nrow - 1, nrow - 1) = (
       C_nonug(nrow - 1, nrow - 1) + s2_nug) * log10;
   }
   if (beta_est) {
     for (int i = 0; i < nrow - 1; i++) {
       for (int j = i + 1; j < nrow; j++) {
         // double total = 1;
         NumericVector dvec(nsum), rvec(nsum);
         for(int k = 0; k < nsum; ++k) {
           // total += theta[k] * pow((x(i,k) - x(j,k)), 2.0);
           double d = fabs(x(i,k) - x(j,k)) / theta[k];
           dvec[k] = d;
           double r = 0;
           if (d <= .5) {
             r = 1-6*pow(d, 2.0)+6*pow(d, 3.0);
           } else if (d <= 1) {
             r = 2*pow(1-d, 3.0);
           } else {
             r = 0;
           }
           rvec[k] = r;
           // total *= r;
         }

         for(int k = 0; k < nsum; ++k) {
           double grad = 0;
           if (x(i,k) - x(j,k) > 0) {
             grad = 1;
           } else {
             grad = -1;
           }
           double d = fabs(x(i,k) - x(j,k)) / theta[k];
           double dr = 0;
           if (d <= .5) {
             // tmp2 = 1-6*pow(d, 2)+6*pow(d,3);
             dr = -12*d+18*pow(d, 2.0);
           } else if (d <= 1) {
             dr = -6*pow(1-d, 2.0);
           } else {
             dr = 0;
           }
           grad *= log10 * (-(x(i,k) - x(j,k))) / theta[k] * dr;
           // if (d > 0) {
           //   grad *= total / rvec[k];
           // } else {
           for (int l=0; l < nsum; l++) {
             if (k != l) {
               grad *= rvec[l];
             }
           }
           grad *= s2;
           // }
           dC_dparams(k,i,j) = grad;
           dC_dparams(k,j,i) = dC_dparams(k,i,j);
         }


       }
     }

     for (int k=0; k < nsum; k++) {
       for (int i = 0; i < nrow; i++) { //# Get diagonal set to zero
         dC_dparams(k,i,i) = 0;
       }
     }
   }

   return dC_dparams;
 }

