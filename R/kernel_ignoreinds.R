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
#' @field kernel Kernel to use on indices that aren't ignored
#' @field ignoreinds Indices to ignore. For a matrix X, these are the columns
#' to ignore. For example, when those dimensions will be given a different
#' kernel, such as for factors.
#' @examples
#' kg <- Gaussian$new(D=3)
#' kig <- GauPro::IgnoreIndsKernel$new(k = Gaussian$new(D=3), ignoreinds = 2)
#' Xtmp <- as.matrix(expand.grid(1:2, 1:2, 1:2))
#' cbind(Xtmp, kig$k(Xtmp))
#' cbind(Xtmp, kg$k(Xtmp))
IgnoreIndsKernel <- R6::R6Class(
  classname = "GauPro_kernel_IgnoreInds",
  inherit = GauPro_kernel,
  active = list(
    #' @field s2_est
    #' Is s2 being estimated?
    s2_est = function(val) {
      if (missing(val)) {
        return(self$kernel$s2_est)
      } else {
        stopifnot(length(val) == 1, is.logical(val))
        if (val) { # Turn on s2_est
          self$kernel$s2_est <- TRUE
        } else { # Turn off s2_est
          self$kernel$s2_est <- FALSE
        }
      }
    },
    #' @field s2
    #' Value of s2 (variance)
    s2 = function(x) {
      self$kernel$s2
    }
  ),
  public = list(
    D = NULL,
    kernel = NULL,
    ignoreinds = NULL,
    #' @description Initialize kernel object
    #' @param k Kernel to use on the non-ignored indices
    #' @param ignoreinds Indices of columns of X to ignore.
    #' @param useC Should C code used? Not implemented for IgnoreInds.
    initialize = function(k, ignoreinds, useC=TRUE) {
      stopifnot("GauPro_kernel" %in% class(k))
      stopifnot(is.numeric(ignoreinds))
      stopifnot(abs(ignoreinds - round(ignoreinds)) < 1e-8, ignoreinds>.99999)
      self$kernel <- k
      self$ignoreinds <- ignoreinds
      self$D <- self$kernel$D + length(self$ignoreinds)
      self$useC <- useC # Never used
    },
    #' @description Calculate covariance between two points
    #' @param x vector.
    #' @param y vector, optional. If excluded, find correlation
    #' of x with itself.
    #' @param ... Passed to kernel
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
    #' @description Find covariance of two points
    #' @param x vector
    #' @param y vector
    #' @param ... Passed to kernel
    kone = function(x, y, ...) {
      x2 = x[-self$ignoreinds]
      y2 = y[-self$ignoreinds]
      self$kernel$k(x=x2, y=y2, ...)
    },
    #' @description Derivative of covariance with respect to parameters
    #' @param params Kernel parameters
    #' @param X matrix of points in rows
    #' @param ... Passed to kernel
    dC_dparams = function(params=NULL, X, ...) {
      X2 <- X[, -self$ignoreinds, drop=FALSE]
      arglist <- list(params=params, X=X2, ...)
      do.call(self$kernel$dC_dparams, arglist)
    },
    #' @description Calculate covariance matrix and its derivative
    #'  with respect to parameters
    #' @param params Kernel parameters
    #' @param X matrix of points in rows
    #' @param nug Value of nugget
    C_dC_dparams = function(params=NULL, X, nug) {
      X2 <- X[, -self$ignoreinds, drop=FALSE]
      arglist <- list(params=params, X=X2, nug=nug)
      do.call(self$kernel$C_dC_dparams, arglist)
    },
    # Below is updated version using arma, was called dC_dx_arma before
    #' @description Derivative of covariance with respect to X
    #' @param XX matrix of points
    #' @param X matrix of points to take derivative with respect to
    #' @param ... Additional arguments passed on to the kernel
    dC_dx = function(XX, X, ...) {
      if (!is.matrix(XX)) {stop("XX must be matrix")}
      if (ncol(X) != ncol(XX)) {stop("XX and X must have same number of cols")}
      # corr_gauss_dCdX(XX, X, theta, s2)
      # stop()
      out <- array(data=0, dim=c(nrow(XX), ncol(X), nrow(X)))
      # useinds <- setdiff(1:D, self$ignoreinds)
      sub_dC_dx <- self$kernel$dC_dx(XX=XX[, -self$ignoreinds, drop=FALSE],
                                     X=X[, -self$ignoreinds, drop=FALSE],
                                     ...)
      out[, -self$ignoreinds, ] <- sub_dC_dx
      out
    },
    #' @description Starting point for parameters for optimization
    #' @param ... Passed to kernel
    param_optim_start = function(...) {
      self$kernel$param_optim_start(...)
    },
    #' @description Starting point for parameters for optimization
    #' @param ... Passed to kernel
    param_optim_start0 = function(...) {
      self$kernel$param_optim_start0(...)
    },
    #' @description Lower bounds of parameters for optimization
    #' @param ... Passed to kernel
    param_optim_lower = function(...) {
      self$kernel$param_optim_lower(...)
    },
    #' @description Upper bounds of parameters for optimization
    #' @param ... Passed to kernel
    param_optim_upper = function(...) {
      self$kernel$param_optim_upper(...)
    },
    #' @description Set parameters from optimization output
    #' @param ... Passed to kernel
    set_params_from_optim = function(...) {
      self$kernel$set_params_from_optim(...)
    },
    #' @description Get s2 from params vector
    #' @param ... Passed to kernel
    s2_from_params = function(...) {
      self$kernel$s2_from_params(...)
    },
    #' @description Print this object
    print = function() {
      p1 <- capture.output(self$kernel)
      p1b <- paste0('\t', p1)
      p1c <- paste0(p1b, collapse = '\n')
      cat('GauPro kernel: Ignore indices\n')
      cat('\tignore inds = ', self$ignoreinds, '\n')
      cat(p1c, "\n")
    }
  ),
  private = list(

  )
)
