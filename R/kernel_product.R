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
#' @field s2 Variance
#' @examples
#' k1 <- Exponential$new(beta=1)
#' k2 <- Matern32$new(beta=2)
#' k <- k1 + k2
#' k$k(matrix(c(2,1), ncol=1))
kernel_product <- R6::R6Class(classname = "GauPro_kernel_product",
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
      self$s2 <- self$k1$s2 * self$k2$s2

      # if (self$k1$s2_est && )
      s2_est <- (self$k1$s2_est || self$k2$s2_est)
    },
    #' @description Calculate covariance between two points
    #' @param x vector.
    #' @param y vector, optional. If excluded, find correlation
    #' of x with itself.
    #' @param params parameters to use instead of beta and s2.
    #' @param ... Not used
    k = function(x, y=NULL, params, ...) {
      if (missing(params)) {
        self$k1$k(x=x, y=y) * self$k2$k(x=x, y=y)
      } else {
        params1 <- params[1:self$k1pl]
        params2 <- params[(self$k1pl+1):(self$k1pl+self$k2pl)]
        self$k1$k(x=x, y=y, params=params1) * self$k2$k(x=x, y=y, params=params2)
      }
    },
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
      s2_1 <- self$k1$s2_from_params(params1)
      s2_2 <- self$k2$s2_from_params(params2)
      # #
      # out1 <- self$k1$dC_dparams(params=params1, C=C, X=X, C_nonug=C_nonug)
      # out2 <- self$k2$dC_dparams(params=params2, C=C, X=X, C_nonug=C_nonug)
      # Can't pass in C, no longer specific to each one
      out1 <- self$k1$C_dC_dparams(params=params1, X=X, nug=s2_2 * nug)
      out2 <- self$k2$C_dC_dparams(params=params2, X=X, nug=s2_1 * nug)
      C1_nonug <- (out1[[1]] - diag(s2_1 * s2_2*nug, nrow(X)))
      C2_nonug <- (out2[[1]] - diag(s2_1 * s2_2*nug, nrow(X)))
      C <- C1_nonug * C2_nonug + diag(s2_1*s2_2*nug, nrow(X))

      # Multiply beta params by opposite C_nonug
      n_beta1 <- length(params1) - self$k1$s2_est
      # browser()
      if (n_beta1 > 0) { # At least 1 beta param
        # out1[[2]][[1:n_beta1]] <- lapply(out1[[2]][1:n_beta1], function(m) {m * C2_nonug})
        for (i in 1:n_beta1) {
          out1[[2]][i,,] <- out1[[2]][i,,] * C2_nonug
        }
      }
      n_beta2 <- length(params2) - self$k2$s2_est
      if (n_beta2 > 0) { # At least 1 beta param
        # out2[[2]][[1:n_beta2]] <- lapply(out2[[2]][1:n_beta2], function(m) {m * C1_nonug})
        for (i in 1:n_beta2) {
          out2[[2]][i,,] <- out2[[2]][i,,] * C1_nonug
        }
      }

      # Fix s2
      if (self$k1$s2_est) { # Fix here, don't like this
        out1[[2]][length(params1),,] <- C * log(10) # on log scale/ s2_1
      }
      if (self$k2$s2_est) { # Fix here, don't like this
        out2[[2]][length(params2),,] <- C * log(10) # on log scale/ s2_2
      }
      dim1 <- dim(out1[[2]])
      dim2 <- dim(out2[[2]])
      dC_dparams <- array(dim=c(dim1[1]+dim2[1], dim1[2], dim1[3]))
      dC_dparams[1:dim1[1],,] <- out1[[2]]
      dC_dparams[(1+dim1[1]):(dim1[1]+dim2[1]),,] <- out2[[2]]

      dC_dparams #c(out1[[2]],out2[[2]]))#, c(out1[[2]]*out2[[2]]))
    },
    #' @description Calculate covariance matrix and its derivative
    #'  with respect to parameters
    #' @param params Kernel parameters
    #' @param X matrix of points in rows
    #' @param nug Value of nugget
    C_dC_dparams = function(params=NULL, X, nug) {#browser(text = "Make sure all in one list")
      params1 <- params[1:self$k1pl]
      params2 <- params[(self$k1pl+1):(self$k1pl+self$k2pl)]
      s2_1 <- self$k1$s2_from_params(params1)
      s2_2 <- self$k2$s2_from_params(params2)
      # #
      # out1 <- self$k1$dC_dparams(params=params1, C=C, X=X, C_nonug=C_nonug)
      # out2 <- self$k2$dC_dparams(params=params2, C=C, X=X, C_nonug=C_nonug)
      # Can't pass in C, no longer specific to each one
      out1 <- self$k1$C_dC_dparams(params=params1, X=X, nug=s2_2 * nug)
      out2 <- self$k2$C_dC_dparams(params=params2, X=X, nug=s2_1 * nug)
      C1_nonug <- (out1[[1]] - diag(s2_1 * s2_2*nug, nrow(X)))
      C2_nonug <- (out2[[1]] - diag(s2_1 * s2_2*nug, nrow(X)))
      C <- C1_nonug * C2_nonug + diag(s2_1*s2_2*nug, nrow(X))

      # Multiply beta params by opposite C_nonug
      n_beta1 <- length(params1) - self$k1$s2_est
      if (n_beta1 > 0) { # At least 1 beta param
        # out1[[2]][[1:n_beta1]] <- lapply(out1[[2]][1:n_beta1], function(m) {m * C2_nonug})
        for (i in 1:n_beta1) {
          out1[[2]][i,,] <- out1[[2]][i,,] * C2_nonug
        }
      }
      n_beta2 <- length(params2) - self$k2$s2_est
      if (n_beta2 > 0) { # At least 1 beta param
        # out2[[2]][[1:n_beta2]] <- lapply(out2[[2]][1:n_beta2], function(m) {m * C1_nonug})
        for (i in 1:n_beta2) {
          out2[[2]][i,,] <- out2[[2]][i,,] * C1_nonug
        }
      }

      # Fix s2
      if (self$k1$s2_est) { # Fix here, don't like this
        out1[[2]][length(params1),,] <- C * log(10) # on log scale/ s2_1
      }
      if (self$k2$s2_est) { # Fix here, don't like this
        out2[[2]][length(params2),,] <- C * log(10) # on log scale/ s2_2
      }
      dim1 <- dim(out1[[2]])
      dim2 <- dim(out2[[2]])
      dC_dparams <- array(dim=c(dim1[1]+dim2[1], dim1[2], dim1[3]))
      dC_dparams[1:dim1[1],,] <- out1[[2]]
      dC_dparams[(1+dim1[1]):(dim1[1]+dim2[1]),,] <- out2[[2]]

      # list(C=C, dC_dparams=c(out1[[2]],out2[[2]]))#, c(out1[[2]]*out2[[2]]))
      list(C=C, dC_dparams=dC_dparams)
    },
    #' @description Derivative of covariance with respect to X
    #' @param XX matrix of points
    #' @param X matrix of points to take derivative with respect to
    dC_dx = function(XX, X) {
      C1 <- self$k1$k(x=XX, y=X)
      C2 <- self$k2$k(x=XX, y=X)
      dC1 <- self$k1$dC_dx(XX=XX, X=X)
      dC2 <- self$k2$dC_dx(XX=XX, X=X)
      for (i in 1:dim(dC2)[2]) {
        dC2[,i,] <- dC2[,i,] * C1
      }
      for (i in 1:dim(dC1)[2]) {
        dC1[,i,] <- dC1[,i,] * C2
      }
      dC2 + dC1
    },
    #' @description Get s2 from params vector
    #' @param params parameter vector
    #' @param s2_est Is s2 being estimated?
    s2_from_params = function(params, s2_est=self$s2_est) {
      params1 <- params[1:self$k1pl]
      params2 <- params[(self$k1pl+1):(self$k1pl+self$k2pl)]
      self$k1$s2_from_params(params=params1) * self$k2$s2_from_params(params=params2)
    }
  )
)
