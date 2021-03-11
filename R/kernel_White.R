#' White noise Kernel R6 class
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
#' @field s2 variance
#' @field logs2 Log of s2
#' @field logs2_lower Lower bound of logs2
#' @field logs2_upper Upper bound of logs2
#' @field s2_est Should s2 be estimated?
#' @examples
#' k1 <- White$new(s2=1e-8)
White <- R6::R6Class(
  classname = "GauPro_kernel_White",
  inherit = GauPro_kernel,
  public = list(
    s2 = NULL, # variance coefficient to scale correlation matrix to covariance
    logs2 = NULL,
    logs2_lower = NULL,
    logs2_upper = NULL,
    s2_est = NULL,
    #' @description Initialize kernel object
    #' @param s2 Initial variance
    #' @param D Number of input dimensions of data
    #' @param s2_lower Lower bound for s2
    #' @param s2_upper Upper bound for s2
    #' @param s2_est Should s2 be estimated?
    initialize = function(s2=1, D, s2_lower=1e-8, s2_upper=1e8, s2_est=TRUE) {
      self$s2 <- s2
      self$logs2 <- log(s2, 10)
      self$logs2_lower <- log(s2_lower, 10)
      self$logs2_upper <- log(s2_upper, 10)
      self$s2_est <- s2_est

      if (missing(D)) {self$D <- NA}
      else {self$D <- D}

    },
    #' @description Calculate covariance between two points
    #' @param x vector.
    #' @param y vector, optional. If excluded, find correlation
    #' of x with itself.
    #' @param s2 Variance parameter.
    #' @param params parameters to use instead of beta and s2.
    k = function(x, y=NULL, s2=self$s2, params=NULL) {#browser()
      if (!is.null(params)) {
        logs2 <- params #[lenpar]
        s2 <- 10^logs2
      } else {
        if (is.null(s2)) {s2 <- self$s2}
      }
      if (is.null(y)) {
        if (is.matrix(x)) {
          val <- diag(s2, nrow(x))
          return(val)
        } else {
          return(s2 * 1)
        }
      }
      if (is.matrix(x) & is.matrix(y)) {
        # s2 * corr_gauss_matrixC(x, y, theta)
        # outer(1:nrow(x), 1:nrow(y), Vectorize(function(i,j){self$kone(x[i,],y[j,],theta=theta, alpha=alpha, s2=s2)}))
        matrix(0, nrow(x), nrow(y))
      } else if (is.matrix(x) & !is.matrix(y)) {
        # s2 * corr_gauss_matrixvecC(x, y, theta)
        # apply(x, 1, function(xx) {self$kone(xx, y, theta=theta, alpha=alpha, s2=s2)})
        matrix(0, nrow(x), 1)
      } else if (is.matrix(y)) {
        # s2 * corr_gauss_matrixvecC(y, x, theta)
        # apply(y, 1, function(yy) {self$kone(yy, x, theta=theta, alpha=alpha, s2=s2)})
        matrix(0, 1, nrow(y))
      } else {
        # self$kone(x, y, theta=theta, alpha=alpha, s2=s2)
        0
      }
    },
    #' @description Find covariance of two points
    #' @param x vector
    #' @param y vector
    #' @param s2 Variance parameter
    kone = function(x, y, s2) {
      # I don't think this should ever be used.
      stop('kernel_white kone should never be used #2398273')
      # Should this return s2 when all(x==y)?
      return(0)
    },
    #' @description Derivative of covariance with respect to parameters
    #' @param params Kernel parameters
    #' @param X matrix of points in rows
    #' @param C_nonug Covariance without nugget added to diagonal
    #' @param C Covariance with nugget
    #' @param nug Value of nugget
    dC_dparams = function(params=NULL, X, C_nonug, C, nug) {#browser(text = "Make sure all in one list")
      n <- nrow(X)
      lenparams <- length(params)
      if (is.null(params)) {
        logs2 <- self$logs2
      } else {
        logs2 <- params
      }
      log10 <- log(10)
      # logs2 <- params[lenparams]
      s2 <- 10 ^ logs2

      if (missing(C_nonug)) { # Assume C missing too, must have nug
        C_nonug <- self$k(x=X, params=params)
        C <- C_nonug + diag(nug*s2, nrow(C_nonug))
      }

      lenparams_D <- as.integer(self$s2_est)
      dC_dparams <- array(dim=c(lenparams, n, n), data=0)
      if (self$s2_est) {
        dC_dparams[lenparams,,] <- C * log10
      }
      return(dC_dparams)
    },
    #' @description Calculate covariance matrix and its derivative
    #'  with respect to parameters
    #' @param params Kernel parameters
    #' @param X matrix of points in rows
    #' @param nug Value of nugget
    C_dC_dparams = function(params=NULL, X, nug) {
      s2 <- self$s2_from_params(params)
      C_nonug <- self$k(x=X, params=params)
      C <- C_nonug + diag(s2*nug, nrow(X))
      dC_dparams <- self$dC_dparams(params=params, X=X, C_nonug=C_nonug, C=C, nug=nug)
      list(C=C, dC_dparams=dC_dparams)
    },
    #' @description Derivative of covariance with respect to X
    #' @param XX matrix of points
    #' @param X matrix of points to take derivative with respect to
    #' @param theta Correlation parameters
    #' @param beta log of theta
    #' @param s2 Variance parameter
    dC_dx = function(XX, X, s2=self$s2) {#browser()
      if (!is.matrix(XX)) {stop()}
      d <- ncol(XX)
      if (ncol(X) != d) {stop()}
      n <- nrow(X)
      nn <- nrow(XX)
      dC_dx <- array(0, dim=c(nn, d, n))
      dC_dx
    },
    #' @description Starting point for parameters for optimization
    #' @param jitter Should there be a jitter?
    #' @param y Output
    #' @param s2_est Is s2 being estimated?
    param_optim_start = function(jitter=F, y, s2_est=self$s2_est) {
      # Use current values for theta, partial MLE for s2
      # vec <- c(log(self$theta, 10), log(sum((y - mu) * solve(R, y - mu)) / n), 10)
      # if (beta_est) {vec <- c(self$beta)} else {vec <- c()}
      # if (alpha_est) {vec <- c(vec, self$logalpha)} else {}
      if (s2_est) {vec <- self$logs2} else {c()}
      # if (jitter && beta_est) {
      #   # vec <- vec + c(self$beta_optim_jitter,  0)
      #   vec[1:length(self$beta)] = vec[1:length(self$beta)] + rnorm(length(self$beta), 0, 1)
      # }
      vec
    },
    #' @description Starting point for parameters for optimization
    #' @param jitter Should there be a jitter?
    #' @param y Output
    #' @param s2_est Is s2 being estimated?
    param_optim_start0 = function(jitter=F, y, s2_est=self$s2_est) {
      # Use 0 for theta, partial MLE for s2
      # vec <- c(rep(0, length(self$theta)), log(sum((y - mu) * solve(R, y - mu)) / n), 10)
      # if (beta_est) {vec <- rep(0, self$beta_length)} else {vec <- c()}
      # if (alpha_est) {vec <- c(vec, 1)} else {}
      if (s2_est) {vec <- 0} else {c()}
      # if (jitter && beta_est) {
      #   vec[1:length(self$beta)] = vec[1:length(self$beta)] + rnorm(length(self$beta), 0, 1)
      # }
      vec
    },
    #' @description Lower bounds of parameters for optimization
    #' @param s2_est Is s2 being estimated?
    param_optim_lower = function(s2_est=self$s2_est) {
      # c(self$beta_lower, self$logs2_lower)
      # if (beta_est) {vec <- c(self$beta_lower)} else {vec <- c()}
      # if (alpha_est) {vec <- c(vec, self$logalpha_lower)} else {}
      if (s2_est) {vec <- self$logs2_lower} else {c()}
      vec
    },
    #' @description Upper bounds of parameters for optimization
    #' @param s2_est Is s2 being estimated?
    param_optim_upper = function(s2_est=self$s2_est) {
      # c(self$beta_upper, self$logs2_upper)
      # if (beta_est) {vec <- c(self$beta_upper)} else {vec <- c()}
      # if (alpha_est) {vec <- c(vec, self$logalpha_upper)} else {}
      if (s2_est) {vec <- self$logs2_upper} else {c()}
      vec
    },
    #' @description Set parameters from optimization output
    #' @param optim_out Output from optimization
    #' @param s2_est s2 estimate
    set_params_from_optim = function(optim_out, s2_est=self$s2_est) {
      loo <- length(optim_out)
      if (s2_est) {
        if (loo != 1) {stop("Error in white kernel #923585")}
        self$logs2 <- optim_out[loo]
        self$s2 <- 10 ^ self$logs2
      } else {
        if (loo != 0) {stop("Error in white kernel #3245235")}
      }
    },
    #' @description Get s2 from params vector
    #' @param params parameter vector
    #' @param s2_est Is s2 being estimated?
    s2_from_params = function(params, s2_est=self$s2_est) {
      # 10 ^ params[length(params)]
      if (s2_est && !is.null(params)) { # Is last if in params
        10 ^ params
      } else { # Else it is just using set value, not being estimated
        self$s2
      }
    }
  )
)
