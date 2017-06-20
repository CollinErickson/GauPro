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



#' Gaussian Kernel R6 class
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
#' k1 <- Gaussian$new(theta=1)
Gaussian <- R6::R6Class(classname = "GauPro_kernel_Gaussian",
  inherit = GauPro_kernel,
  public = list(
    theta = NULL,
    theta_lower = NULL,
    theta_upper = NULL,
    s2 = NULL, # variance coefficient to scale correlation matrix to covariance
    initialize = function(theta, s2=1, theta_lower=0, theta_upper=1e6) {
      self$theta <- theta
      # if (length(theta) == 1) {
      #   self$theta <- rep(theta, self$d)
      # }
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
    l = function(X, y, theta, s2, mu, n) {
      R <- self$r(X, theta)
      n*log(s2) + log(det(R)) + sum(y - mu, Rinv %*% (y-mu))
    },
    dl_dthetas2 = function(X, y, theta, mu, s2, n, firstiter) {
      R <- self$r(X, theta)
      dl_ds2 <- n / s2 - s2^2 * sum((y - mu) * solve(R, y - mu))
      # p should be theta length
      dl_dt <- sapply(1:self$p, function(l) {
        # dR_dti <- R
        dr_dtl <- outer(1:n, 1:n, function(i, j) {-(X[i,k] - X[j,k])^2 * R[i,j]})
        dR_dtl_Rinv <- solve(dR_dtl, R)
        dl_dtl <- diag(dR_dtl) / s2 + sum(Rinv %*% (y-mu), dR_dtl %*% (y-mu))/ s2^2
        dl_dtl
      })
      c(cl_dtl, dl_ds2)
    },
    optim_param_start = function(random, y) {
      if (random) {
        c(log(self$theta, 10) + rnorm(self$p, 0, 1), log(sum((y - mu) * solve(R, y - mu)) / n), 10)
      } else {
        c(log(self$theta, 10), log(sum((y - mu) * solve(R, y - mu)) / n), 10)
      }
    },
    optim_param_lower = function() {
      c(rep(-5, self$p), 8)
    },
    optim_param_upper = function() {
      c(rep(5, self$p), 8)
    },
    optim_fngr = function(X, y, params, mu, n) {
      theta <- 10^params[1:self$p]
      s2 <- 10^params[self$p+1]
      list(fn=self$l(X=X, y=y, theta=theta, s2=s2, mu=mu, n=n),
           gr=self$dl_dthetas2(X=X, y=y, theta=theta, s2=s2, mu=mu, n=n, firstiter=FALSE)
      )
    },
    param_set = function(optim_out) {
      self$theta <- 10^optim_out[1:self$p]
      self$s2 <- 10^optim_out[self$p+1]
    }
  ),
  private = list(

  )
)
