# Trend functions should implement:
# mu prediction for new matrix
# update_params
# get_optim_functions: return optim.func, optim.grad, optim.fngr
# param_optim_lower - lower bound of params
# param_optim_upper - upper
# param_optim_start - current param values
# param_optim_start0 - some central param values that can be used for optimization restarts
# param_optim_jitter - how to jitter params in optimization




#' Trend R6 class
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
#' t1 <- trend_LM$new(D=2)
trend_LM <- R6::R6Class(
  classname = "GauPro_trend_LM",
  inherit = GauPro_trend,
  public = list(
    m = NULL,
    m_lower = NULL,
    m_upper = NULL,
    m_est = NULL,
    b = NULL,
    b_lower = NULL,
    b_upper = NULL,
    b_est = NULL,
    D = NULL,
    initialize = function(D, m = rep(0,D), m_lower=rep(-Inf,D), m_upper=rep(Inf,D), m_est=rep(TRUE,D),
                          b = 0, b_lower=-Inf, b_upper=Inf, b_est=TRUE) {
      self$D <- D
      self$m <- m
      self$m_lower <- m_lower
      self$m_upper <- m_upper
      self$m_est <- m_est
      self$b <- b
      self$b_lower <- b_lower
      self$b_upper <- b_upper
      self$b_est <- b_est
    },
    Z = function(X, m=self$m, b=self$b, params=NULL) {#browser()
      if (!is.null(params)) {m <- params[2:(self$D+1)]; b <- params[1]}
      if (is.matrix(X)) {
        b + X %*% m
        # matrix(m, nrow=nrow(X), ncol=1)
      } else { # If a vector then just a single value
        b + sum(X * m)
      }
    },
    dZ_dparams = function(X, m=self$m_est, b=self$b_est, params=NULL) {
      # Gradient is -2 * t(yminusmu) %*% Siginv %*% du/db
      if (!is.null(params)) {m <- params[2:(self$D+1)]; b <- params[1]}
      # matrix(1, nrow=nrow(X), ncol=1)
      cbind(rep(1,nrow(X)), X)
    },
    dZ_dx = function(X, m=self$m, params=NULL) {
      if (!is.null(params)) {m <- params}
      # matrix(0, nrow=nrow(X), ncol=1)
      matrix(m, nrow=nrow(X), ncol=length(m), byrow=T)

    },
    param_optim_start = function(jitter, trend_est) {
      c(self$b, self$m)
    },
    param_optim_start0 = function(jitter, trend_est) {
      c(self$b, self$m)
    },
    param_optim_lower = function(jitter, trend_est) {
      c(self$b_lower, self$m_lower)
    },
    param_optim_upper = function(jitter, trend_est) {
      c(self$b_upper, self$m_upper)
    },
    set_params_from_optim = function(optim_out) {
      self$b <- optim_out[1]
      self$m <- optim_out[2:(self$D+1)]
    }
  )
)
