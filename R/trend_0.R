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
# @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @field m Trend parameters
#' @field m_lower m lower bound
#' @field m_upper m upper bound
#' @field m_est Should m be estimated?
#' @examples
#' t1 <- trend_0$new()
trend_0 <- R6::R6Class(classname = "GauPro_trend_0",
  inherit = GauPro_trend,
  public = list(
    m = NULL,
    m_lower = NULL,
    m_upper = NULL,
    m_est = NULL,
    #' @description Initialize trend object
    #' @param D Number of input dimensions of data
    #' @param m trend initial parameters
    #' @param m_lower trend lower bounds
    #' @param m_upper trend upper bounds
    #' @param m_est Logical of whether each param should be estimated
    initialize = function(m = 0, m_lower=0, m_upper=0, m_est=FALSE, D=NA) {
      self$m <- m
      self$m_lower <- m_lower
      self$m_upper <- m_upper
      self$m_est <- m_est
      self$D <- D
    },
    #' @description Get trend value for given matrix X
    #' @param X matrix of points
    #' @param m trend parameters
    #' @param params trend parameters
    Z = function(X, m=self$m, params=NULL) {
      if (!is.null(params)) {m <- params}
      if (is.matrix(X)) {
        # matrix(m, nrow=nrow(X), ncol=1)
        matrix(0, nrow(X), 1)
      } else { # If a vector then just a single value
        # m
        c()
      }
    },
    #' @description Derivative of trend with respect to trend parameters
    #' @param X matrix of points
    #' @param m trend values
    #' @param params overrides m
    dZ_dparams = function(X, m=m$est, params=NULL) {
      # Gradient is -2 * t(yminusmu) %*% Siginv %*% du/db
      if (!is.null(params)) {m <- params}
      matrix(0, nrow=nrow(X), ncol=0)
      array(0, dim = c(nrow(X), ncol=0))
    },
    #' @description Derivative of trend with respect to X
    #' @param X matrix of points
    #' @param m trend values
    #' @param params overrides m
    dZ_dx = function(X, m=self$m, params=NULL) {
      if (!is.null(params)) {m <- params}
      matrix(0, nrow=nrow(X), ncol=1)

    },
    #' @description Get parameter initial point for optimization
    #' @param jitter Not used
    #' @param trend_est If the trend should be estimate.
    param_optim_start = function(jitter, trend_est) {
      c()
    },
    #' @description Get parameter initial point for optimization
    #' @param jitter Not used
    #' @param trend_est If the trend should be estimate.
    param_optim_start0 = function(jitter, trend_est) {
      c()
    },
    #' @description Get parameter lower bounds for optimization
    #' @param jitter Not used
    #' @param trend_est If the trend should be estimate.
    param_optim_lower = function(jitter, trend_est) {
      c()
    },
    #' @description Get parameter upper bounds for optimization
    #' @param jitter Not used
    #' @param trend_est If the trend should be estimate.
    param_optim_upper = function(jitter, trend_est) {
      c()
    },
    #' @description Set parameters after optimization
    #' @param optim_out Output from optim
    set_params_from_optim = function(optim_out) {
      # Nothing to set
    }
  )
)
