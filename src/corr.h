using namespace Rcpp;
using namespace arma;

#ifndef PKG_FOO1_H
#define PKG_FOO1_H

Rcpp::NumericMatrix corr_gauss_matrix_symC(NumericMatrix x, NumericVector theta);

arma::mat corr_gauss_matrix_sym_armaC(arma::mat x, arma::vec theta);

#endif
