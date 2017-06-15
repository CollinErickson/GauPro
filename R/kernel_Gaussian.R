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
    initialize = function(theta) {
      self$theta <- theta
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

    }
  ),
  private = list(

  )
)
