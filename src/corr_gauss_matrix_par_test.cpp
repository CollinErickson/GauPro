#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

struct corr_gauss_matrix_par_struct : public Worker {

  // input matrix to read from
  const RMatrix<double> mat1;
  const RMatrix<double> mat2;
  const RVector<double> theta;

  // output matrix to write to
  RMatrix<double> rmat;

  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  corr_gauss_matrix_par_struct(const NumericMatrix mat1,
                               const NumericMatrix mat2,
                               const NumericVector theta,
                               NumericMatrix rmat)
    : mat1(mat1), mat2(mat2), theta(theta), rmat(rmat) {}

  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    int ncol = mat2.ncol();
    double tsum;
    for (std::size_t i = begin; i < end; i++) {
      for (std::size_t j = 0; j < mat2.nrow(); j++) {

        // rows we will operate on
        // RMatrix<double>::Row row1 = mat.row(i);
        // RMatrix<double>::Row row2 = mat.row(j);

        // compute the average using std::tranform from the STL
        // std::vector<double> avg(row1.length());
        // std::transform(row1.begin(), row1.end(), // input range 1
        //                row2.begin(),             // input range 2
        //                avg.begin(),              // output range
        //                average);                 // function to apply

        // calculate divergences
        // double d1 = kl_divergence(row1.begin(), row1.end(), avg.begin());
        // double d2 = kl_divergence(row2.begin(), row2.end(), avg.begin());

        // write to output matrix
        // rmat(i,j) = sqrt(.5 * (d1 + d2));
        // rmat(i,j) = exp(-sum(theta * pow(mat1.row(i) - mat2.row(j), 2)));
        // rmat(i,j) = exp(-sum(theta * pow(mat1.row(i) - mat2.row(j), 2)));
        tsum=0;
        for (int k=0; k < ncol; k++) {
          tsum += theta[k] * pow(mat1(i,k) - mat2(j,k), 2);
        }
        rmat(i,j) = exp(- tsum);
      }
    }
  }
};

//' Correlation Gaussian matrix in C using RcppParallel
//'
//' Faster than nonparallel version for D < 12 and > 20 rows
//' @param x Matrix x
//' @param y Matrix y, must have same number of columns as x
//' @param theta Theta vector
//' @return Correlation matrix
//' @examples
//' corr_gauss_matrixC(matrix(c(1,0,0,1),2,2), matrix(c(1,0,1,1),2,2), c(1,1))
//' @export
// [[Rcpp::export]]
NumericMatrix corr_gauss_matrixCpar(NumericMatrix mat1, NumericMatrix mat2,
                                    NumericVector theta) {

  // allocate the matrix we will return
  NumericMatrix rmat(mat1.nrow(), mat2.nrow());

  // create the worker
  corr_gauss_matrix_par_struct corr_gauss_matrix_par_instance(mat1, mat2,
                                                              theta, rmat);

  // call it with parallelFor
  parallelFor(0, mat1.nrow(), corr_gauss_matrix_par_instance, 20);

  return rmat;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
*/
