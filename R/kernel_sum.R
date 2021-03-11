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
#' @field k1 kernel 1
#' @field k2 kernel 2
#' @field k1_param_length param length of kernel 1
#' @field k2_param_length param length of kernel 2
#' @field k1pl param length of kernel 1
#' @field k2pl param length of kernel 2
#' @field s2 variance
#' @examples
#' k1 <- Exponential$new(beta=1)
#' k2 <- Matern32$new(beta=2)
#' k <- k1 + k2
#' k$k(matrix(c(2,1), ncol=1))
kernel_sum <- R6::R6Class(classname = "GauPro_kernel_sum",
  inherit = GauPro_kernel,
  public = list(
    k1 = NULL,
    k2 = NULL,
    k1_param_length = NULL,
    k2_param_length = NULL,
    k1pl = NULL,
    k2pl = NULL,
    s2 = NULL,
    #' @description Initialize kernel
    #' @param k1 Kernel 1
    #' @param k2 Kernel 2
    initialize = function(k1, k2) {
      self$k1 <- k1
      self$k2 <- k2
      self$k1_param_length <- length(self$k1$param_optim_start())
      self$k1pl <- self$k1_param_length
      self$k2_param_length <- length(self$k2$param_optim_start())
      self$k2pl <- self$k2_param_length
      self$s2 <- self$k1$s2 + self$k2$s2
    },
    #' @description Calculate covariance between two points
    #' @param x vector.
    #' @param y vector, optional. If excluded, find correlation
    #' of x with itself.
    #' @param params parameters to use instead of beta and s2.
    #' @param ... Not used
    k = function(x, y=NULL, params, ...) {
      if (missing(params)) {
        self$k1$k(x=x, y=y) + self$k2$k(x=x, y=y)
      } else {
        params1 <- params[1:self$k1pl]
        params2 <- params[(self$k1pl+1):(self$k1pl+self$k2pl)]
        self$k1$k(x=x, y=y, params=params1) + self$k2$k(x=x, y=y, params=params2)
      }
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
    #' @description Starting point for parameters for optimization
    #' @param jitter Should there be a jitter?
    #' @param y Output
    param_optim_start = function(jitter=F, y) {
      # Use current values for theta, partial MLE for s2
      # vec <- c(log(self$theta, 10), log(sum((y - mu) * solve(R, y - mu)) / n), 10)
      c(self$k1$param_optim_start(jitter=jitter), self$k2$param_optim_start(jitter=jitter))
    },
    #' @description Starting point for parameters for optimization
    #' @param jitter Should there be a jitter?
    #' @param y Output
    param_optim_start0 = function(jitter=F, y) {
      # Use 0 for theta, partial MLE for s2
      # vec <- c(rep(0, length(self$theta)), log(sum((y - mu) * solve(R, y - mu)) / n), 10)
      c(self$k1$param_optim_start0(jitter=jitter), self$k2$param_optim_start0(jitter=jitter))
    },
    #' @description Lower bounds of parameters for optimization
    param_optim_lower = function() {
      c(self$k1$param_optim_lower(), self$k2$param_optim_lower())
    },
    #' @description Upper bounds of parameters for optimization
    param_optim_upper = function() {
      c(self$k1$param_optim_upper(), self$k2$param_optim_upper())
    },
    #' @description Set parameters from optimization output
    #' @param optim_out Output from optimization
    set_params_from_optim = function(optim_out) {
      oo1 <- optim_out[1:self$k1pl]
      self$k1$set_params_from_optim(optim_out=oo1)
      oo2 <- optim_out[(self$k1pl+1):(self$k1pl+self$k2pl)]
      self$k2$set_params_from_optim(optim_out=oo2)
      self$s2 <- self$k1$s2 + self$k2$s2
    },
    #' @description Derivative of covariance with respect to parameters
    #' @param params Kernel parameters
    #' @param X matrix of points in rows
    #' @param C_nonug Covariance without nugget added to diagonal
    #' @param C Covariance with nugget
    #' @param nug Value of nugget
    dC_dparams = function(params=NULL, C, X, C_nonug, nug) {#browser(text = "Make sure all in one list")
      params1 <- params[1:self$k1pl]
      params2 <- params[(self$k1pl+1):(self$k1pl+self$k2pl)]
      # #
      # out1 <- self$k1$dC_dparams(params=params1, C=C, X=X, C_nonug=C_nonug)
      # out2 <- self$k2$dC_dparams(params=params2, C=C, X=X, C_nonug=C_nonug)
      # Can't pass in C, no longer specific to each one
      out1 <- self$k1$dC_dparams(params=params1, X=X, nug=nug)
      out2 <- self$k2$dC_dparams(params=params2, X=X, nug=nug)
      dim1 <- dim(out1)
      dim2 <- dim(out2)
      dC_dparams <- array(dim=c(dim1[1]+dim2[1], dim1[2], dim1[3]))
      dC_dparams[1:dim1[1],,] <- out1
      dC_dparams[(1+dim1[1]):(dim1[1]+dim2[1]),,] <- out2
      dC_dparams
      # list(c(out1[[1]],out2[[1]]), c(out1[[2]]+out2[[2]]))
    },
    #' @description Calculate covariance matrix and its derivative
    #'  with respect to parameters
    #' @param params Kernel parameters
    #' @param X matrix of points in rows
    #' @param nug Value of nugget
    C_dC_dparams = function(params=NULL, X, nug) {#browser(text = "Make sure all in one list")
      params1 <- params[1:self$k1pl]
      params2 <- params[(self$k1pl+1):(self$k1pl+self$k2pl)]
      # out1 <- self$k1$dC_dparams(params=params1, C=C, X=X, C_nonug=C_nonug)
      # out2 <- self$k2$dC_dparams(params=params2, C=C, X=X, C_nonug=C_nonug)
      # list(c(out1[[1]],out2[[1]]), c(out1[[2]]+out2[[2]]))
      # cat('In kernel_sum C_dC_params\n')

      # Need to recalculate C for each so pass nug instead
      out1 <- self$k1$C_dC_dparams(params=params1, X=X, nug=nug)
      out2 <- self$k2$C_dC_dparams(params=params2, X=X, nug=nug)
      C <- out1[[1]] + out2[[1]]
      # dC_dparams <- c(out1[[2]],out2[[2]])#, c(out1[[2]]+out2[[2]])
      dim1 <- dim(out1[[2]])
      dim2 <- dim(out2[[2]])
      dC_dparams <- array(dim=c(dim1[1]+dim2[1], dim1[2], dim1[3]))
      dC_dparams[1:dim1[1],,] <- out1[[2]]
      dC_dparams[(1+dim1[1]):(dim1[1]+dim2[1]),,] <- out2[[2]]
      # dC_dparams
      list(C=C, dC_dparams=dC_dparams)
    },
    #' @description Derivative of covariance with respect to X
    #' @param XX matrix of points
    #' @param X matrix of points to take derivative with respect to
    dC_dx = function(XX, X) {
      self$k1$dC_dx(XX=XX, X=X) + self$k2$dC_dx(XX=XX, X=X)
    },
    #' @description Get s2 from params vector
    #' @param params parameter vector
    #' @param s2_est Is s2 being estimated?
    s2_from_params = function(params) {
      params1 <- params[1:self$k1pl]
      params2 <- params[(self$k1pl+1):(self$k1pl+self$k2pl)]
      self$k1$s2_from_params(params=params1) + self$k2$s2_from_params(params=params2)
    }
  )
)
