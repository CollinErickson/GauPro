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

// This was slower, vectorized version below is same speed or 15% slower in 2D,
//   3x faster in 6D, 6.6x faster in 12D
// //' Correlation Gaussian matrix in C using Rcpp
// //' @param x Matrix x
// //' @param y Matrix y, must have same number of columns as x
// //' @param theta Theta vector
// //' @return Correlation matrix
// //' @examples
// //' corr_gauss_matrixC(matrix(c(1,0,0,1),2,2), matrix(c(1,0,1,1),2,2), c(1,1))
// //' @export
// // [[Rcpp::export]]
// NumericMatrix corr_gauss_matrixC(NumericMatrix x, NumericMatrix y, NumericVector theta) {
//   int nrow = x.nrow(), ncol = y.nrow();
//   int nsum = x.ncol();
//   NumericMatrix out(nrow, ncol);
//
//   for (int i = 0; i < nrow; i++) {
//     for (int j = 0; j < ncol; j++) {
//
//       double total = 0;
//       for(int k = 0; k < nsum; ++k) {
//         total += theta[k] * pow((x(i,k) - y(j,k)), 2);
//       }
//       total = exp(-total);
//
//       out(i, j) = total;
//     }
//   }
//   return out;
// }

//' Correlation Gaussian matrix in C using Rcpp
//' @param x Matrix x
//' @param y Matrix y, must have same number of columns as x
//' @param theta Theta vector
//' @return Correlation matrix
//' @examples
//' corr_gauss_matrixC(matrix(c(1,0,0,1),2,2), matrix(c(1,0,1,1),2,2), c(1,1))
//' @export
// [[Rcpp::export]]
NumericMatrix corr_gauss_matrixC(NumericMatrix x, NumericMatrix y, NumericVector theta) {
  int nrow = x.nrow(), ncol = y.nrow();
  // int ndim = x.ncol();
  NumericMatrix out(nrow, ncol);
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      // Vectorized way with sugar, should work
      out(i, j) = std::exp(-sum(theta * Rcpp::pow(x.row(i) - y.row(j), 2.0)));
      // // As a for loop.
      // out(i, j) = 0;
      // for (int k=0; k < ndim; k++) {
      //   out(i,j) -= theta[k] * std::pow(x(i, k) - y(j, k), 2);
      // }
      // out(i, j) = std::exp(out(i, j));
    }
  }
  return out;
}

// Tried to make faster for 2D, failed
// //' Correlation Gaussian matrix in C using Rcpp
// //' @param x Matrix x
// //' @param y Matrix y, must have same number of columns as x
// //' @param theta Theta vector
// //' @return Correlation matrix
// //' @examples
// //' corr_gauss_matrixC(matrix(c(1,0,0,1),2,2), matrix(c(1,0,1,1),2,2), c(1,1))
// //' @export
// // [[Rcpp::export]]
// NumericMatrix corr_gauss_matrixC2D(NumericMatrix x, NumericMatrix y, NumericVector theta) {
//   int nrow = x.nrow(), ncol = y.nrow();
//   NumericMatrix out(nrow, ncol);
//   for (int i = 0; i < nrow; i++) {
//     for (int j = 0; j < ncol; j++) {
//       // out(i, j) = exp(-sum(theta * pow(x.row(i) - y.row(j), 2.0)));
//       // out(i, j) = exp(- theta(0)*pow(x(i,0)-y(j,0), 2.0) - theta(1) * pow(x(i,1) - y(j,1), 2.0));
//       out(i, j) = exp(- theta(0)*(x(i,0)-y(j,0))*(x(i,0)-y(j,0)) - theta(1) * (x(i,1)-y(j,1))*(x(i,1)-y(j,1)));
//     }
//   }
//   return out;
// }


//' Correlation Gaussian matrix in C (symmetric)
//' @param x Matrix x
//' @param theta Theta vector
//' @return Correlation matrix
//' @export
//' @examples
//' corr_gauss_matrix_symC(matrix(c(1,0,0,1),2,2),c(1,1))
// [[Rcpp::export]]
NumericMatrix corr_gauss_matrix_symC(NumericMatrix x, NumericVector theta) {
  int nrow = x.nrow();
  int nsum = x.ncol();
  NumericMatrix out(nrow, nrow);

  for (int i = 0; i < nrow - 1; i++) {
    for (int j = i + 1; j < nrow; j++) {

      double total = 0;
      for(int k = 0; k < nsum; ++k) {
        total += theta[k] * std::pow((x(i,k) - x(j,k)), 2);
      }
      total = std::exp(-total);

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
NumericVector corr_gauss_matrixvecC(NumericMatrix x, NumericVector y, NumericVector theta) {
  int nrow = x.nrow(); //, ncol = y.nrow();
  int nsum = x.ncol();
  NumericVector out(nrow);

  for (int i = 0; i < nrow; i++) {
    double total = 0;
    for(int k = 0; k < nsum; ++k) {
      total += theta[k] * std::pow((x(i,k) - y(k)), 2);
    }
    total = std::exp(-total);

    out(i) = total;
  }
  return out;
}



//' Correlation Gaussian matrix in C using Armadillo (symmetric)
//'
//' About 30% faster than Rcpp version.
//' @param x Matrix x
//' @param theta Theta vector
//' @return Correlation matrix
//' @examples
//' corr_gauss_matrix_sym_armaC(matrix(c(1,0,0,1),2,2),c(1,1))
//'
//' x3 <- matrix(runif(1e3*6), ncol=6)
//' th <- runif(6)
//' t3 <- corr_gauss_matrix_symC(x3, th)
//' t4 <- corr_gauss_matrix_sym_armaC(x3, th)
//' identical(t3, t4)
//' # microbenchmark::microbenchmark(corr_gauss_matrix_symC(x3, th),
//' #                     corr_gauss_matrix_sym_armaC(x3, th), times=50)
//' @export
// [[Rcpp::export]]
arma::mat corr_gauss_matrix_sym_armaC(arma::mat x, arma::vec theta) {
  int nrow = x.n_rows;
  int nsum = x.n_cols;
  arma::mat out(nrow, nrow);

  for (int i = 0; i < nrow - 1; i++) {
    for (int j = i + 1; j < nrow; j++) {

      double total = 0;
      for(int k = 0; k < nsum; ++k) {
        total += theta[k] * std::pow((x(i,k) - x(j,k)), 2);
      }
      total = std::exp(-total);

      out(i, j) = total;
      out(j, i) = total;
    }
  }
  for (int i = 0; i < nrow; i++) {
    out(i, i) = 1;
  }
  return out;
}


//' Correlation Gaussian matrix in C using Armadillo
//'
//' 20-25% faster than Rcpp version.
//' @param x Matrix x
//' @param y Matrix y, must have same number of columns as x
//' @param theta Theta vector
//' @param s2 Variance to multiply matrix by
//' @return Correlation matrix
//' @examples
//' corr_gauss_matrix_armaC(matrix(c(1,0,0,1),2,2),matrix(c(1,0,1,1),2,2),c(1,1))
//'
//' x1 <- matrix(runif(100*6), nrow=100, ncol=6)
//' x2 <- matrix(runif(1e4*6), ncol=6)
//' th <- runif(6)
//' t1 <- corr_gauss_matrixC(x1, x2, th)
//' t2 <- corr_gauss_matrix_armaC(x1, x2, th)
//' identical(t1, t2)
//' # microbenchmark::microbenchmark(corr_gauss_matrixC(x1, x2, th),
//' #                                corr_gauss_matrix_armaC(x1, x2, th))
//' @export
// [[Rcpp::export]]
arma::mat corr_gauss_matrix_armaC(arma::mat x, arma::mat y, arma::vec theta, double s2 = 1.0) {
  int nrowx = x.n_rows;
  int nrowy = y.n_rows;
  int nsum = x.n_cols;
  arma::mat out(nrowx, nrowy, arma::fill::zeros);
  for(int k = 0; k < nsum; ++k) {
    for (int i = 0; i < nrowx; i++) {
      for (int j = 0; j < nrowy; j++) {
        out(i,j) += theta[k] * std::pow((x(i,k) - y(j,k)), 2);
      }
    }
  }
  out = arma::exp(-out);
  if (s2 != 1.0) {
    out *= s2;
  }
  return out;
}




// Trying to get C_dC for Gaussian kernel
//' Correlation Gaussian matrix in C (symmetric)
//' @param x Matrix x
//' @param theta Theta vector
//' @param C_nonug cov mat without nugget
//' @param s2_est whether s2 is being estimated
//' @param beta_est Whether theta/beta is being estimated
//' @param lenparams_D Number of parameters the derivative is being calculated for
//' @param s2_nug s2 times the nug
//' @return Correlation matrix
//' @export
//' @examples
//' corr_gauss_matrix_symC(matrix(c(1,0,0,1),2,2),c(1,1))
// [[Rcpp::export]]
arma::cube kernel_gauss_dC(arma::mat x, arma::vec theta, arma::mat C_nonug, bool s2_est, bool beta_est, int lenparams_D, double s2_nug) {
  int nrow = x.n_rows;
  int nsum = x.n_cols;
  arma::cube dC_dparams(lenparams_D, nrow, nrow);

  //arma::mat out(nrow, nrow)
  /*for (int i = 0; i < nrow - 1; i++) {
    for (int j = i + 1; j < nrow; j++) {

      double total = 0;
      for(int k = 0; k < nsum; ++k) {
        total += theta[k] * pow((x(i,k) - x(j,k)), 2);
      }
      total = exp(-total);

      out(i, j) = total;
      out(j, i) = total; // since symmetric
    }
  }
  for (int i = 0; i < nrow; i++) {
    out(i, i) = 1;
  }*/

  if (s2_est) {
    // dC_dparams(lenparams_D,,) = C * log(10.0);
    for (int i = 0; i < nrow - 1; i++) {
      for (int j = i + 1; j < nrow; j++) {
        dC_dparams(lenparams_D - 1,i,j) = C_nonug(i,j) * log(10.0);
        dC_dparams(lenparams_D - 1,j,i) = dC_dparams(lenparams_D - 1,i,j);
      }
      dC_dparams(lenparams_D - 1, i, i) = (C_nonug(i,i) + s2_nug) * log(10.0);
    }
    dC_dparams(lenparams_D - 1, nrow - 1, nrow - 1) = (C_nonug(nrow - 1, nrow - 1) + s2_nug) * log(10.0);
  }
  if (beta_est) {
    for (int k = 0; k < nsum; k++) {
      for (int i = 0; i < nrow - 1; i++) {
        for (int j = i + 1; j < nrow; j++) {
          dC_dparams(k,i,j) = - C_nonug(i,j) * std::pow(x(i,k) - x(j,k), 2) * theta(k) * log(10.0);
          dC_dparams(k,j,i) = dC_dparams(k,i,j);
        }
      }

      for (int i = 0; i < nrow; i++) { //# Get diagonal set to zero
        dC_dparams(k,i,i) = 0;
      }
    }
  }

  return dC_dparams;
}





// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
#timesTwo(42)
*/
