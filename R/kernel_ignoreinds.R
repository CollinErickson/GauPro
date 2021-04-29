#' Kernel R6 class
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
#' @field D Number of input dimensions of data
#' @examples
#' kg <- Gaussian$new(D=3)
#' kig <- GauPro::IgnoreIndsKernel$new(k = Gaussian$new(D=3), ignoreinds = 2)
#' Xtmp <- as.matrix(expand.grid(1:2, 1:2, 1:2))
#' cbind(Xtmp, kig$k(Xtmp))
#' cbind(Xtmp, kg$k(Xtmp))
IgnoreIndsKernel <- R6::R6Class(
  classname = "GauPro_kernel_IgnoreInds",
  active = list(
    s2_est = function(x) {
      self$kernel$s2_est
    },
    s2 = function(x) {
      self$kernel$s2
    }
  ),
  public = list(
    D = NULL,
    kernel = NULL,
    ignoreinds = NULL,
    initialize = function(k, ignoreinds) {
      stopifnot("GauPro_kernel" %in% class(k))
      self$kernel <- k
      self$ignoreinds <- ignoreinds
    },
    k = function(x, y=NULL, ...) {
      if (is.matrix(x)) {
        x2 = x[, -self$ignoreinds, drop=FALSE]
      } else {
        x2 = x[-self$ignoreinds]
      }
      if (!is.null(y)) {
        if (is.matrix(y)) {
          y2 = y[, -self$ignoreinds, drop=FALSE]
        } else {

          y2 <- y[-self$ignoreinds]
        }
      } else {
        y2 <- NULL
      }
      self$kernel$k(x=x2, y=y2, ...)
    },
    kone = function(x, y, ...) {
      x2 = x[-self$ignoreinds]
      y2 = y[-self$ignoreinds]
      self$kernel$k(x=x2, y=y2, ...)
    },
    dC_dparams = function(params=NULL, X, ...) {
      X2 <- X[, -self$ignoreinds, drop=FALSE]
      arglist <- list(params=params, X=X2, ...)
      do.call(self$kernel$dC_dparams, arglist)
    },
    C_dC_dparams = function(params=NULL, X, nug) {
      X2 <- X[, -self$ignoreinds, drop=FALSE]
      arglist <- list(params=params, X=X2, nug=nug)
      do.call(self$kernel$C_dC_dparams, arglist)
    },
    param_optim_start = function(...) {
      self$kernel$param_optim_start(...)
    },
    param_optim_start0 = function(...) {
      self$kernel$param_optim_start0(...)
    },
    param_optim_lower = function(...) {
      self$kernel$param_optim_lower(...)
    },
    param_optim_upper = function(...) {
      self$kernel$param_optim_upper(...)
    },
    set_params_from_optim = function(...) {
      self$kernel$set_params_from_optim(...)
    },
    s2_from_params = function(...) {
      self$kernel$s2_from_params(...)
    }

  ),
  private = list(

  )
)
