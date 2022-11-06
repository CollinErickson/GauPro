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
# @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' k1 <- Gaussian_l$new(l=1)
Gaussian_l <- R6::R6Class(classname = "GauPro_kernel_Gaussian_lengthscale",
  inherit = GauPro_kernel,
  public = list(
    l = NULL,
    # l_lower = NULL,
    # l_upper = NULL,
    logl_lower = NULL,
    logl_upper = NULL,
    l_length = NULL,
    s2 = NULL, # variance coefficient to scale correlation matrix to covariance
    logs2 = NULL,
    logs2_lower = NULL,
    logs2_upper = NULL,
    initialize = function(l, s2=1, l_lower=1e-8, l_upper=1e6,
                          s2_lower=1e-8, s2_upper=1e8) {
      if (any(l <= 0)) {stop("l must be > 0")}
      self$l <- l
      self$l_length <- length(l)
      # if (length(theta) == 1) {
      #   self$theta <- rep(theta, self$d)
      # }
      self$logl_lower <- log(l_lower, 10)
      self$logl_upper <- log(l_upper, 10)

      self$s2 <- s2
      self$logs2 <- log(s2, 10)
      self$logs2_lower <- log(s2_lower, 10)
      self$logs2_upper <- log(s2_upper, 10)
    },
    k = function(x, y=NULL, l=self$l, s2=self$s2, params=NULL) {
      if (!is.null(params)) {
        logl <- params[1:(length(params)-1)]
        l <- 10 ^ logl
        logs2 <- params[length(params)]
        s2 <- 10 ^ logs2
      } else {
        if (is.null(l)) {l <- self$l}
        if (is.null(s2)) {s2 <- self$s2}
      }
      theta <- .5 / l^2
      if (is.null(y)) {
        if (is.matrix(x)) {
          cgmtry <- try(val <- s2 * corr_gauss_matrix_symC(x, theta))
          if (inherits(cgmtry,"try-error")) {stop("Error cgmtry")}
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
    # k1 = function(x, y, l=self$l) {
    #   theta <- .5 / l^2
    #   self$s2 * exp(-sum(theta * (x-y)^2))
    # },
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
      vec <- c(log(self$l, 10), self$logs2)
      if (jitter) {
        # vec <- vec + c(self$beta_optim_jitter,  0)
        vec[1:length(self$l)] = vec[1:length(self$l)] + rnorm(length(self$l), 0, 1)
      }
      vec
    },
    param_optim_start0 = function(jitter=F, y) {
      # Use 0 for theta, partial MLE for s2
      # vec <- c(rep(0, length(self$theta)), log(sum((y - mu) * solve(R, y - mu)) / n), 10)
      vec <- c(rep(0, self$l_length), 0)
      if (jitter) {
        vec[1:length(self$l)] = vec[1:length(self$l)] + rnorm(length(self$l), 0, 1)
      }
      vec
    },
    param_optim_lower = function() {
      c(self$logl_lower, self$logs2_lower)
    },
    param_optim_upper = function() {
      c(self$logl_upper, self$logs2_upper)
    },
    set_params_from_optim = function(optim_out) {
      loo <- length(optim_out)
      self$l <- 10 ^ optim_out[1:(loo-1)]
      self$logs2 <- optim_out[loo]
      self$s2 <- 10 ^ self$logs2
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
    dC_dparams = function(params=NULL, C, X, C_nonug, nug) {#browser(text = "Make sure all in one list")
      if (is.null(params)) {params <- c(log(self$l,10), self$logs2)}
      lenpar <- length(params)
      logl <- params[1:(lenpar - 1)]
      l <- 10 ^ logl
      one_over_l2 <- 1 / l^2
      theta <- .5 * one_over_l2 #/ l^2
      log10 <- log(10)
      logs2 <- params[lenpar]
      s2 <- 10 ^ logs2
      dC_dlogs2 <- C * log10 # / s2
      dC_dls <- rep(list(C_nonug), length(l))
      n <- nrow(X)
      for (k in 1:length(l)) {
        for (i in seq(1, n-1, 1)) {
          for (j in seq(i+1, n, 1)) {
            dC_dls[[k]][i,j] <- dC_dls[[k]][i,j] * (X[i,k] - X[j,k])^2 * one_over_l2[k] * log10
            dC_dls[[k]][j,i] <- dC_dls[[k]][i,j]
          }
        }
        for (i in seq(1, n, 1)) { # Get diagonal set to zero
          dC_dls[[k]][i,i] <- 0
        }
      }

      mats <- c(dC_dls, list(dC_dlogs2))
      return(list(dC_dparams=mats,
                  s2
      ))
    },
    # param_set = function(optim_out) {
    #   # self$theta <- 10^optim_out[1:self$p]
    #   # self$s2 <- 10^optim_out[self$p+1]
    #   self$l <- optim_out[1:self$l_length]
    #   self$s2 <- optim_out[self$l_length+1]
    # },
    s2_from_params = function(params) {
      10 ^ params[length(params)]
    }
  ),
  private = list(

  )
)
