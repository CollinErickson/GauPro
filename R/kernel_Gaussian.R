# Kernels should implement:
# k kernel function for two vectors
# update_params
# get_optim_functions: return optim.func, optim.grad, optim.fngr
# param_optim_lower - lower bound of params
# param_optim_upper - upper
# param_optim_start - current param values
# param_optim_start0 - some central param values that can be used for optimization restarts
# param_optim_jitter - how to jitter params in optimization

# Suggested
# deviance
# deviance_grad
# deviance_fngr
# grad



#' Kernel R6 class
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @useDynLib GauPro, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats optim
#' @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#'
Gaussian <- R6::R6Class(classname = "GauPro_kernel_Gaussian",
  inherit = GauPro_kernel,
  public = list(
    theta = NULL,
    s2 = NULL, # variance coefficient to scale correlation matrix to covariance
    initialize = function(theta, s2=1, theta_lower=0, theta_upper=1e6) {
      self$theta <- theta
      if (length(theta) == 1) {
        self$theta <- rep(theta, self$d)
      }
      self$s2 <- s2
      self$theta_lower <- theta_lower
      self$theta_upper <- theta_upper
    },
    k = function(x, y=NULL, theta=self$theta) {
      if (is.null(y)) {
        if (is.matrix(x)) {
          return(self$s2 * corr_gauss_matrix_symC(x, theta))
        } else {
          return(self$s2 * 1)
        }
      }
      if (is.matrix(x) & is.matrix(y)) {
        self$s2 * corr_gauss_matrixC(x, y, theta)
      } else if (is.matrix(x) & !is.matrix(y)) {
        self$s2 * corr_gauss_matrixvecC(x, y, theta)
      } else if (is.matrix(y)) {
        self$s2 * corr_gauss_matrixvecC(y, x, theta)
      } else {
        self$s2 * exp(-sum(theta * (x-y)^2))
      }
    },
    k1 = function(x, y, theta=self$theta) {
      self$s2 * exp(-sum(theta * (x-y)^2))
    },
    km = function(x, y, theta=self$theta) {

    },
    dl_dthetas2 = function(X, y, theta, mu, s2, n, firstiter) {
      R <- self$r(X, theta)
      dl_ds2 <- n / s2 - s2^2 * sum((y - u) * solve(R, y - mu))
      # p should be theta length
      dl_dt <- sapply(1:self$p, function(l) {
        # dR_dti <- R
        dr_dtl <- outer(1:n, 1:n, function(i, j) {-(X[i,k] - X[j,k])^2 * R[i,j]})
        dR_dtl_Rinv <- solve(dR_dtl, R)
        dl_dtl <- diag() / s2 + sum(Rinv %*% (y-u), dR_dtl %*% (y-u))/ s2^2
        dl_dtl
      })
      c(cl_dtl, dl_ds2)
    }
  ),
  private = list(

  )
)
