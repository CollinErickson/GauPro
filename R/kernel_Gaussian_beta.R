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
#' k1 <- Gaussian_beta$new(beta=1)
Gaussian_beta <- R6::R6Class(classname = "GauPro_kernel_Gaussian_beta",
  inherit = GauPro_kernel,
  public = list(
    beta = NULL,
    beta_lower = NULL,
    beta_upper = NULL,
    beta_length = NULL,
    s2 = NULL, # variance coefficient to scale correlation matrix to covariance
    s2_lower = NULL,
    s2_upper = NULL,
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
      self$s2_lower <- s2_lower
      self$s2_upper <- s2_upper
    },
    k = function(x, y=NULL, beta=self$beta, s2=self$s2, params=NULL) {
      if (!is.null(params)) {
        beta <- params[1:(length(params)-1)]
        s2 <- params[length(params)]
      } else {#browser()
        if (is.null(beta)) {theta <- self$beta}
        if (is.null(s2)) {s2 <- self$s2}
      }
      theta <- 10^beta
      if (is.null(y)) {
        if (is.matrix(x)) {#browser()
          cgmtry <- try(val <- s2 * corr_gauss_matrix_symC(x, theta))
          if (inherits(cgmtry,"try-error")) {browser()}
          return(val)
        } else {
          return(s2 * 1)
        }
      }
      if (is.matrix(x) & is.matrix(y)) {
        s2 * corr_gauss_matrixC(x, y, theta)
      } else if (is.matrix(x) & !is.matrix(y)) {
        s2 * corr_gauss_matrixvecC(x, y, theta)
      } else if (is.matrix(y)) {
        s2 * corr_gauss_matrixvecC(y, x, theta)
      } else {
        s2 * exp(-sum(theta * (x-y)^2))
      }
    },
    k1 = function(x, y, beta=self$beta) {
      theta <- 10^beta
      self$s2 * exp(-sum(theta * (x-y)^2))
    },
    # l = function(X, y, beta, s2, mu, n) {
    #   theta <- 10^beta
    #   R <- self$r(X, theta)
    #   n*log(s2) + log(det(R)) + sum(y - mu, Rinv %*% (y-mu))
    # },
    # dl_dbetas2 = function(X, y, beta, mu, s2, n, firstiter) {
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
    # beta_optim_jitter = function() {
    #   rnorm(self$p, 0, 1)
    # },
    param_optim_start = function(jitter=F, y) {
      # Use current values for theta, partial MLE for s2
      # vec <- c(log(self$theta, 10), log(sum((y - mu) * solve(R, y - mu)) / n), 10)
      vec <- c(self$beta, self$s2)
      if (jitter) {
        # vec <- vec + c(self$beta_optim_jitter,  0)
        vec[1:length(self$beta)] = vec[1:length(self$beta)] + rnorm(length(self$beta), 0, 1)
      }
      vec
    },
    param_optim_start0 = function(jitter=F, y) {
      # Use 0 for theta, partial MLE for s2
      # vec <- c(rep(0, length(self$theta)), log(sum((y - mu) * solve(R, y - mu)) / n), 10)
      vec <- c(rep(1, self$beta_length), 1)
      if (jitter) {
        vec[1:length(self$beta)] = vec[1:length(self$beta)] + rnorm(length(self$beta), 0, 1)
      }
      vec
    },
    param_optim_lower = function() {
      c(self$beta_lower, self$s2_lower)
    },
    param_optim_upper = function() {
      c(self$beta_upper, self$s2_upper)
    },
    set_params_from_optim = function(optim_out) {
      loo <- length(optim_out)
      self$beta <- optim_out[1:(loo-1)]
      self$s2 <- optim_out[loo]
    },
    # optim_fngr = function(X, y, params, mu, n) {
    #   theta <- 10^params[1:self$p]
    #   s2 <- 10^params[self$p+1]
    #   list(fn=self$l(X=X, y=y, theta=theta, s2=s2, mu=mu, n=n),
    #        gr=self$dl_dthetas2(X=X, y=y, theta=theta, s2=s2, mu=mu, n=n, firstiter=FALSE)
    #   )
    # },
    # get_optim_functions = function(param_update) {
    #
    # },
    dC_dparams = function(params=NULL, C, X, C_nonug) {#browser(text = "Make sure all in one list")
      if (is.null(params)) {params <- c(self$beta, self$s2)}
      beta <- params[1:(length(params) - 1)]
      theta <- 10^beta
      log10 <- log(10)
      s2 <- tail(params, 1)
      dC_ds2 <- C / s2
      dC_dbetas <- rep(list(C_nonug), length(beta))
      n <- nrow(X)
      for (k in 1:length(beta)) {
        for (i in seq(1, n-1, 1)) {
          for (j in seq(i+1, n, 1)) {
            dC_dbetas[[k]][i,j] <- - dC_dbetas[[k]][i,j] * (X[i,k] - X[j,k])^2 * theta[k] * log10
            dC_dbetas[[k]][j,i] <- dC_dbetas[[k]][i,j]
          }
        }
        for (i in seq(1, n, 1)) { # Get diagonal set to zero
          dC_dbetas[[k]][i,i] <- 0
        }
      }

      mats <- c(dC_dbetas, list(dC_ds2))
      return(list(dC_dparams=mats,
                  s2
      ))
    },
    param_set = function(optim_out) {
      # self$theta <- 10^optim_out[1:self$p]
      # self$s2 <- 10^optim_out[self$p+1]
      self$beta <- optim_out[1:self$beta_length]
      self$s2 <- optim_out[self$beta_length+1]
    },
    s2_from_params = function(params=params) {
      params[length(params)]
    }
  ),
  private = list(

  )
)
