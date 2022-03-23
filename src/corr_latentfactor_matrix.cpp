//#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
//using namespace arma;


//' Correlation Latent factor  matrix in C (symmetric)
//' @param x Matrix x
//' @param theta Theta vector
//' @param xindex Index to use
//' @param latentdim Number of latent dimensions
//' @param offdiagequal What to set off-diagonal values with matching values to.
//' @return Correlation matrix
//' @export
//' @examples
//' corr_latentfactor_matrix_symC(matrix(c(1,.5, 2,1.6, 1,0),ncol=2,byrow=T), c(1.5,1.8), 1, 1, 1-1e-6)
//' corr_latentfactor_matrix_symC(matrix(c(0,0,0,1,0,0,0,2,0,0,0,3,0,0,0,4), ncol=4, byrow=T),
//'   c(0.101, -0.714, 0.114, -0.755, 0.117, -0.76, 0.116, -0.752),
//'   4, 2, 1-1e-6) * 6.85
// [[Rcpp::export]]
NumericMatrix corr_latentfactor_matrix_symC(NumericMatrix x, NumericVector theta,
                                            int xindex, int latentdim,
                                            double offdiagequal) {
  int nrow = x.nrow();
  int nsum = x.ncol();
  NumericMatrix out(nrow, nrow);
  int xindoffset, yindoffset;
  int xlev;
  int ylev;
  double total;
  for (int i = 0; i < nrow - 1; i++) {
    for (int j = i + 1; j < nrow; j++) {
      //Rcout << "i "  << i << " j " << j << "\n";
      xlev = x(i, xindex - 1);
      ylev = x(j, xindex - 1);
      //Rcout << "xlev "  << xlev << " ylev " << ylev << "\n";
      if (xlev == ylev) {
        total = offdiagequal;
      } else {
        xindoffset = (xlev - 1) * latentdim;
        yindoffset = (ylev - 1) * latentdim;
        //Rcout << "offsets are " << xindoffset << yindoffset << i << j << "\n";
        //Rcout << "i vals " << i << xindex << latentdim << xindoffset << "\n";
        //Rcout << " j vals" << j << xindex << latentdim << yindoffset << "\n";
        //Rcout << "x matrix:" <<  x << "\n";
        total = 0;
        double latx, laty;
        for(int k = 0; k < latentdim; ++k) {
          //total += theta[k] * pow((x(i,k) - x(j,k)), 2.0);
          latx = theta[xindoffset + k];
          laty = theta[yindoffset + k];
          total += pow(latx - laty, 2);
        }
        total = exp(-total);
      }
      //total = sqrt(3 * total);
      //total = (1 + total) * exp(-total);
      //total = x[i,j];
      out(i, j) = total;
      out(j, i) = total; // since symmetric
    }
  }
  for (int i = 0; i < nrow; i++) {
    out(i, i) = 1;
  }
  return out;
}

/*** R
corr_latentfactor_matrix_symC(matrix(c(0,0,0,1,0,0,0,2,0,0,0,3,0,0,0,4), ncol=4, byrow=T),
                                c(0.101, -0.714, 0.114, -0.755, 0.117, -0.76, 0.116, -0.752),
                                4, 2, 1-1e-6) * 6.85
*/

