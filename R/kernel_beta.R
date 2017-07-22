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



#' Beta Kernel R6 class
#'
#' This is the base structure for a kernel that uses beta = log10(theta)
#' for the lengthscale parameter.
#' It standardizes the params because they all use the same underlying structure.
#' Kernels that inherit this only need to implement kone and dC_dparams.
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
#' k1 <- Matern52$new(beta=0)
GauPro_kernel_beta <- R6::R6Class(classname = "GauPro_kernel_beta",
  inherit = GauPro_kernel,
  public = list(
    beta = NULL,
    beta_lower = NULL,
    beta_upper = NULL,
    beta_length = NULL,
    s2 = NULL, # variance coefficient to scale correlation matrix to covariance
    logs2 = NULL,
    logs2_lower = NULL,
    logs2_upper = NULL,
    initialize = function(beta, s2=1, beta_lower=-8, beta_upper=6,
                          s2_lower=1e-8, s2_upper=1e8) {
      self$beta <- beta
      self$beta_length <- length(beta)
      # if (length(theta) == 1) {
      #   self$theta <- rep(theta, self$d)
      # }
      self$beta_lower <- beta_lower
      self$beta_upper <- beta_upper

      self$s2 <- s2
      self$logs2 <- log(s2, 10)
      self$logs2_lower <- log(s2_lower, 10)
      self$logs2_upper <- log(s2_upper, 10)
    },
    k = function(x, y=NULL, beta=self$beta, s2=self$s2, params=NULL) {#browser()
      if (!is.null(params)) {
        lenpar <- length(params)
        beta <- params[1:(lenpar-1)]
        logs2 <- params[lenpar]
        s2 <- 10^logs2
      } else {#browser()
        if (is.null(beta)) {beta <- self$beta}
        if (is.null(s2)) {s2 <- self$s2}
      }
      theta <- 10^beta
      if (is.null(y)) {
        if (is.matrix(x)) {#browser()
          # cgmtry <- try(val <- s2 * corr_gauss_matrix_symC(x, theta))
          val <- outer(1:nrow(x), 1:nrow(x), Vectorize(function(i,j){self$kone(x[i,],x[j,],theta=theta, s2=s2)}))
          # if (inherits(cgmtry,"try-error")) {browser()}
          return(val)
        } else {
          return(s2 * 1)
        }
      }
      if (is.matrix(x) & is.matrix(y)) {
        # s2 * corr_gauss_matrixC(x, y, theta)
        outer(1:nrow(x), 1:nrow(y), Vectorize(function(i,j){self$kone(x[i,],y[j,],theta=theta, s2=s2)}))
      } else if (is.matrix(x) & !is.matrix(y)) {
        # s2 * corr_gauss_matrixvecC(x, y, theta)
        apply(x, 1, function(xx) {self$kone(xx, y, theta=theta, s2=s2)})
      } else if (is.matrix(y)) {
        # s2 * corr_gauss_matrixvecC(y, x, theta)
        apply(y, 1, function(yy) {self$kone(yy, x, theta=theta, s2=s2)})
      } else {
        self$kone(x, y, theta=theta, s2=s2)
      }
    },
    kone = function(x, y, beta, theta, s2) {
      # Kernels that inherit should implement this or k.
    },
    param_optim_start = function(jitter=F, y) {
      # Use current values for theta, partial MLE for s2
      # vec <- c(log(self$theta, 10), log(sum((y - mu) * solve(R, y - mu)) / n), 10)
      vec <- c(self$beta, self$logs2)
      if (jitter) {
        # vec <- vec + c(self$beta_optim_jitter,  0)
        vec[1:length(self$beta)] = vec[1:length(self$beta)] + rnorm(length(self$beta), 0, 1)
      }
      vec
    },
    param_optim_start0 = function(jitter=F, y) {
      # Use 0 for theta, partial MLE for s2
      # vec <- c(rep(0, length(self$theta)), log(sum((y - mu) * solve(R, y - mu)) / n), 10)
      vec <- c(rep(0, self$beta_length), 0)
      if (jitter) {
        vec[1:length(self$beta)] = vec[1:length(self$beta)] + rnorm(length(self$beta), 0, 1)
      }
      vec
    },
    param_optim_lower = function() {
      c(self$beta_lower, self$logs2_lower)
    },
    param_optim_upper = function() {
      c(self$beta_upper, self$logs2_upper)
    },
    set_params_from_optim = function(optim_out) {
      loo <- length(optim_out)
      self$beta <- optim_out[1:(loo-1)]
      self$logs2 <- optim_out[loo]
      self$s2 <- 10 ^ self$logs2
    },
    # dC_dparams = function(params=NULL, C, X, C_nonug) {
    #   Kernels that inherit from this must implement this.
    # },
    s2_from_params = function(params) {
      10 ^ params[length(params)]
    }
  ),
  private = list(

  )
)
