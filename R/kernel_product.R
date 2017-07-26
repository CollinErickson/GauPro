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
    k = function(x, y=NULL, params, ...) {
      if (missing(params)) {
        self$k1$k(x=x, y=y) + self$k2$k(x=x, y=y)
      } else {
        params1 <- params[1:self$k1pl]
        params2 <- params[(self$k1pl+1):(self$k1pl+self$k2pl)]
        self$k1$k(x=x, y=y, params=params1) * self$k2$k(x=x, y=y, params=params2)
      }
    },
    param_optim_start = function(jitter=F, y) {
      # Use current values for theta, partial MLE for s2
      # vec <- c(log(self$theta, 10), log(sum((y - mu) * solve(R, y - mu)) / n), 10)
      c(self$k1$param_optim_start(jitter=jitter), self$k2$param_optim_start(jitter=jitter))
    },
    param_optim_start0 = function(jitter=F, y) {
      # Use 0 for theta, partial MLE for s2
      # vec <- c(rep(0, length(self$theta)), log(sum((y - mu) * solve(R, y - mu)) / n), 10)
      c(self$k1$param_optim_start0(jitter=jitter), self$k2$param_optim_start0(jitter=jitter))
    },
    param_optim_lower = function() {
      c(self$k1$param_optim_lower(), self$k2$param_optim_lower())
    },
    param_optim_upper = function() {
      c(self$k1$param_optim_upper(), self$k2$param_optim_upper())
    },
    set_params_from_optim = function(optim_out) {
      oo1 <- optim_out[1:self$k1pl]
      self$k1$set_params_from_optim(optim_out=oo1)
      oo2 <- optim_out[(self$k1pl+1):(self$k1pl+self$k2pl)]
      self$k2$set_params_from_optim(optim_out=oo2)
      self$s2 <- self$k1$s2 + self$k2$s2
    },
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
      if (n_beta1 > 0) { # At least 1 beta param
        # out1[[2]][[1:n_beta1]] <- lapply(out1[[2]][1:n_beta1], function(m) {m * C2_nonug})
        for (i in 1:n_beta1) {
          out1[[2]][[i]] <- out1[[2]][[i]] * C2_nonug
        }
      }
      n_beta2 <- length(params2) - self$k2$s2_est
      if (n_beta2 > 0) { # At least 1 beta param
        # out2[[2]][[1:n_beta2]] <- lapply(out2[[2]][1:n_beta2], function(m) {m * C1_nonug})
        for (i in 1:n_beta2) {
          out2[[2]][[i]] <- out2[[2]][[i]] * C1_nonug
        }
      }

      # Fix s2
      if (self$k1$s2_est) { # Fix here, don't like this
        out1[[2]][[length(params1)]] <- C * log(10) # on log scale/ s2_1
      }
      if (self$k2$s2_est) { # Fix here, don't like this
        out2[[2]][[length(params2)]] <- C * log(10) # on log scale/ s2_2
      }

      list(dC_dparams=c(out1[[2]],out2[[2]]))#, c(out1[[2]]*out2[[2]]))
    },
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
          out1[[2]][[i]] <- out1[[2]][[i]] * C2_nonug
        }
      }
      n_beta2 <- length(params2) - self$k2$s2_est
      if (n_beta2 > 0) { # At least 1 beta param
        # out2[[2]][[1:n_beta2]] <- lapply(out2[[2]][1:n_beta2], function(m) {m * C1_nonug})
        for (i in 1:n_beta2) {
          out2[[2]][[i]] <- out2[[2]][[i]] * C1_nonug
        }
      }

      # Fix s2
      if (self$k1$s2_est) { # Fix here, don't like this
        out1[[2]][[length(params1)]] <- C * log(10) # on log scale/ s2_1
      }
      if (self$k2$s2_est) { # Fix here, don't like this
        out2[[2]][[length(params2)]] <- C * log(10) # on log scale/ s2_2
      }

      list(C=C, dC_dparams=c(out1[[2]],out2[[2]]))#, c(out1[[2]]*out2[[2]]))
    },
    s2_from_params = function(params, s2_est=self$s2_est) {
      params1 <- params[1:self$k1pl]
      params2 <- params[(self$k1pl+1):(self$k1pl+self$k2pl)]
      self$k1$s2_from_params(params=params1) * self$k2$s2_from_params(params=params2)
    }
  )
)
