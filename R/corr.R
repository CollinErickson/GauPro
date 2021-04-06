corr_gauss_noC <- function(a, b, theta) {#browser()
  exp(-sum(theta * (a-b)^2))
}

corr_gauss_matrix_noC <- function(x, x2=x, theta) {#browser()
  #outer(x,x2, gauss_cor)
  outer(1:nrow(x),1:nrow(x2), Vectorize(function(i,j) corr_gauss_noC(x[i,], x2[j,], theta=theta)))
}

#' Gaussian correlation
#'
#' @param x First data matrix
#' @param x2 Second data matrix
#' @param theta Correlation parameter
#'
#' @return Correlation matrix
#' @export
#'
#' @examples
#' corr_gauss_matrix(matrix(1:10,ncol=1), matrix(6:15,ncol=1), 1e-2)
corr_gauss_matrix <- function(x, x2=NULL, theta) {
  stopifnot(is.matrix(x), is.vector(theta))
  if (is.null(x2)) {
    corr_gauss_matrix_symC(x, theta)
  } else {
    stopifnot(is.matrix(x2), ncol(x)==ncol(x2), ncol(x)==length(theta))
    corr_gauss_matrixC(x, x2, theta)
  }
}


# corr_gauss using C++
# Rcpp::cppFunction('double corr_gaussC_wrongplace(NumericVector a, NumericVector b, NumericVector theta) {
#   int n = a.size();
#   double total = 0;
#   for(int i = 0; i < n; ++i) {
#   total += theta[i] * pow((a[i] - b[i]), 2.0);
#   }
#   total = exp(-total);
#   return total;
#   }')
if (F) {
  corr_gaussC(1:5, 6:10, 1e-2/(1:5))

  system.time(replicate(1e5, corr_gaussC(1:10, 6:15, 1e-2/(1:10))))
}

# Rcpp::cppFunction('NumericMatrix corr_gauss_matrixC_wrongplace(NumericMatrix x, NumericMatrix y, NumericVector theta) {
#   int nrow = x.nrow(), ncol = y.nrow();
#   int nsum = x.ncol();
#   NumericMatrix out(nrow, ncol);
#
#   for (int i = 0; i < nrow; i++) {
#     for (int j = 0; j < ncol; j++) {
#
#       double total = 0;
#       for(int k = 0; k < nsum; ++k) {
#         total += theta[k] * pow((x(i,k) - y(j,k)), 2.0);
#       }
#       total = exp(-total);
#
#       out(i, j) = total;
#     }
#   }
#   return out;
#   }')
if (F) {
  corr_gauss_matrixC(matrix(1:6,2,3), matrix(1:9,3,3), theta=1:3)
  system.time(replicate(1e5, corr_gauss_matrix(matrix(1:6,2,3), matrix(1:9,3,3), theta=1:3)))
}
