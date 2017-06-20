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
#' k2 <- Gaussian$new(theta=2)
#' k <- k1 + k2
#' k$k(matrix(c(2,1), ncol=1))
kernel_product <- R6::R6Class(classname = "GauPro_kernel_product",
  inherit = GauPro_kernel,
  public = list(
    k1 = NULL,
    k2 = NULL,
    initialize = function(k1, k2) {
      self$k1 <- k1
      self$k2 <- k2
    },
    k = function(x, y=NULL, ...) {
      self$k1$k(x=x, y=y) * self$k2$k(x=x, y=y)
    },
    # k1 = function(x, y, theta=self$theta) {
    #   self$s2 * exp(-sum(theta * (x-y)^2))
    # },
    # km = function(x, y, theta=self$theta) {
    #
    # },
    # l = function(X, y, theta, s2, mu, n) {
    #   R <- self$r(X, theta)
    #   n*log(s2) + log(det(R)) + sum(y - mu, Rinv %*% (y-mu))
    # },
    # dl_dthetas2 = function(X, y, theta, mu, s2, n, firstiter) {
    #   R <- self$r(X, theta)
    #   dl_ds2 <- n / s2 - s2^2 * sum((y - mu) * solve(R, y - mu))
    #   # p should be theta length
    #   dl_dt <- sapply(1:self$p, function(l) {
    #     # dR_dti <- R
    #     dr_dtl <- outer(1:n, 1:n, function(i, j) {-(X[i,k] - X[j,k])^2 * R[i,j]})
    #     dR_dtl_Rinv <- solve(dR_dtl, R)
    #     dl_dtl <- diag(dR_dtl) / s2 + sum(Rinv %*% (y-mu), dR_dtl %*% (y-mu))/ s2^2
    #     dl_dtl
    #   })
    #   c(cl_dtl, dl_ds2)
    # },
    optim_param_start = function(random, y) {
      c(self$k1$optim_param_start(random=random, y=y),
        self$k2$optim_param_start(random=random, y=y))
    },
    optim_param_lower = function() {
      c(self$k1$optim_param_lower(),
        self$k2$optim_param_lower())
    },
    optim_param_upper = function() {
      c(self$k1$optim_param_upper(),
        self$k2$optim_param_upper())
    }
    # optim_fngr = function(X, y, params, mu, n) {
    #   theta <- 10^params[1:self$p]
    #   s2 <- 10^params[self$p+1]
    #   list(fn=self$l(X=X, y=y, theta=theta, s2=s2, mu=mu, n=n),
    #        gr=self$dl_dthetas2(X=X, y=y, theta=theta, s2=s2, mu=mu, n=n, firstiter=FALSE)
    #   )
    # },
    # param_set = function(optim_out) {
    #   self$theta <- 10^optim_out[1:self$p]
    #   self$s2 <- 10^optim_out[self$p+1]
    # }
  )
)
