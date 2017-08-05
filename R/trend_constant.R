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
#' t1 <- trend_c$new()
trend_c <- R6::R6Class(
  classname = "GauPro_trend_c",
  inherit = GauPro_trend,
  public = list(
    m = NULL,
    m_lower = NULL,
    m_upper = NULL,
    m_est = NULL,
    initialize = function(m = 0, m_lower=-Inf, m_upper=Inf, m_est=TRUE) {
      self$m <- m
      self$m_lower <- m_lower
      self$m_upper <- m_upper
      self$m_est <- m_est
    },
    Z = function(X, m=self$m, params=NULL) {
      if (!is.null(params)) {m <- params}
      if (is.matrix(X)) {
        matrix(m, nrow=nrow(X), ncol=1)
      } else { # If a vector then just a single value
        m
      }
    },
    dZ_dparams = function(X, m=self$m, params=NULL) {
      # Gradient is -2 * t(yminusmu) %*% Siginv %*% du/db
      if (!is.null(params)) {m <- params}
      matrix(1, nrow=nrow(X), ncol=1)
      # array(1, dim=c(1, nrow(X), 1))
    },
    dZ_dx = function(X, m=self$m, params=NULL) {
      if (!is.null(params)) {m <- params}
      matrix(0, nrow=nrow(X), ncol=1)

    },
    param_optim_start = function(jitter, trend_est) {
      0
    },
    param_optim_start0 = function(jitter, trend_est) {
      0
    },
    param_optim_lower = function(jitter, trend_est) {
      -Inf
    },
    param_optim_upper = function(jitter, trend_est) {
      Inf
    },
    set_params_from_optim = function(optim_out) {
      self$m <- optim_out
    }
  )
)
