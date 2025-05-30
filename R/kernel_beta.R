#' Beta Kernel R6 class
#'
#' This is the base structure for a kernel that uses beta = log10(theta)
#' for the lengthscale parameter.
#' It standardizes the params because they all use the same underlying
#' structure.
#' Kernels that inherit this only need to implement kone and dC_dparams.
#'
#' @docType class
#' @importFrom R6 R6Class
# @export
#' @useDynLib GauPro, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats optim
# @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link[R6]{R6Class}} with methods for fitting GP model.
#' @format \code{\link[R6]{R6Class}} object.
#' @field beta Parameter for correlation. Log of theta.
#' @field beta_est Should beta be estimated?
#' @field beta_lower Lower bound of beta
#' @field beta_upper Upper bound of beta
#' @field beta_length length of beta
#' @field s2 variance
#' @field logs2 Log of s2
#' @field logs2_lower Lower bound of logs2
#' @field logs2_upper Upper bound of logs2
#' @field s2_est Should s2 be estimated?
#' @field useC Should C code used? Much faster.
#' @field isotropic If isotropic then a single beta/theta is used for all
#' dimensions. If not (anisotropic) then a separate beta/beta is used for
#' each dimension.
#' @examples
#' #k1 <- Matern52$new(beta=0)
GauPro_kernel_beta <- R6::R6Class(
  classname = "GauPro_kernel_beta",
  inherit = GauPro_kernel,
  public = list(
    beta = NULL,
    beta_lower = NULL,
    beta_upper = NULL,
    beta_length = NULL,
    beta_est = NULL, # Should beta be estimated?
    s2 = NULL, # variance coefficient to scale correlation matrix to covariance
    logs2 = NULL,
    logs2_lower = NULL,
    logs2_upper = NULL,
    s2_est = NULL, # Should s2 be estimated?
    useC = NULL,
    isotropic = NULL,
    #' @description Initialize kernel object
    #' @param beta Initial beta value
    #' @param s2 Initial variance
    #' @param D Number of input dimensions of data
    #' @param beta_lower Lower bound for beta
    #' @param beta_upper Upper bound for beta
    #' @param beta_est Should beta be estimated?
    #' @param s2_lower Lower bound for s2
    #' @param s2_upper Upper bound for s2
    #' @param s2_est Should s2 be estimated?
    #' @param useC Should C code used? Much faster.
    #' @param isotropic If isotropic then a single beta/theta is used for all
    #' dimensions. If not (anisotropic) then a separate beta/beta is used for
    #' each dimension.
    initialize = function(beta, s2=1, D,
                          beta_lower=-8, beta_upper=6, beta_est=TRUE,
                          s2_lower=1e-8, s2_upper=1e8, s2_est=TRUE,
                          useC=TRUE, isotropic=FALSE
    ) {
      # Check beta and D
      missing_beta <- missing(beta)
      missing_D    <- missing(D)
      if (missing_beta && missing_D) {stop("Must give kernel beta or D")}
      else if (missing_beta) {beta <- rep(0, D)}
      else if (missing_D) {D <- length(beta)}
      else {if (length(beta) != D) {stop("beta and D should have same length")}}

      self$D <- D
      self$beta <- beta
      self$beta_length <- length(beta)

      # Setting beta_lower so dimensions are right
      self$beta_lower <- if (length(beta_lower) == self$beta_length) {
        beta_lower
      } else if (length(beta_lower)==1) {
        rep(beta_lower,
            self$beta_length)
      } else {
        stop("Error for kernel_beta beta_lower")
      }
      #self$beta_upper <- beta_upper
      self$beta_upper <- if (length(beta_upper) == self$beta_length) {
        beta_upper
      } else if (length(beta_upper)==1) {
        rep(beta_upper, self$beta_length)
      } else {
        stop("Error for kernel_beta beta_upper")
      }
      stopifnot(self$beta_lower <= self$beta_upper)

      stopifnot(is.logical(isotropic), length(isotropic) == 1)
      self$isotropic <- isotropic
      if (isotropic) {
        self$beta <- self$beta[1]
        self$beta_lower <- self$beta_lower[1]
        self$beta_upper <- self$beta_upper[1]
        self$beta_length <- 1
      }

      stopifnot(length(self$beta) == self$beta_length,
                length(self$beta_lower) == self$beta_length,
                length(self$beta_upper) == self$beta_length)

      self$beta_est <- beta_est

      self$s2 <- s2
      self$logs2 <- log(s2, 10)
      stopifnot(s2_lower <= s2, s2 <= s2_upper, s2_lower <= s2_upper)
      self$logs2_lower <- log(s2_lower, 10)
      self$logs2_upper <- log(s2_upper, 10)
      self$s2_est <- s2_est

      self$useC <- useC
    },
    #' @description Calculate covariance between two points
    #' @param x vector.
    #' @param y vector, optional. If excluded, find correlation
    #' of x with itself.
    #' @param beta Correlation parameters. Log of theta.
    #' @param s2 Variance parameter.
    #' @param params parameters to use instead of beta and s2.
    k = function(x, y=NULL, beta=self$beta, s2=self$s2, params=NULL) {
      if (!is.null(params)) {
        lenpar <- length(params)
        beta <- params[1:(lenpar-1)]
        logs2 <- params[lenpar]
        s2 <- 10^logs2
      } else {
        if (is.null(beta)) {beta <- self$beta}
        if (is.null(s2)) {s2 <- self$s2}
      }
      theta <- 10^beta
      if (self$isotropic && length(theta) == self$beta_length) {
        theta <- rep(theta, self$D)
      }
      if (is.null(y)) {
        if (is.matrix(x)) {
          # cgmtry <- try(val <- s2 * corr_gauss_matrix_symC(x, theta))
          val <- outer(1:nrow(x), 1:nrow(x),
                       Vectorize(function(i,j){self$kone(x[i,],x[j,],
                                                         theta=theta, s2=s2)}))
          # if (inherits(cgmtry,"try-error")) {browser()}
          return(val)
        } else {
          return(s2 * 1)
        }
      }
      if (is.matrix(x) & is.matrix(y)) {
        # s2 * corr_gauss_matrixC(x, y, theta)
        outer(1:nrow(x), 1:nrow(y),
              Vectorize(function(i,j){self$kone(x[i,],y[j,],
                                                theta=theta, s2=s2)}))
      } else if (is.matrix(x) & !is.matrix(y)) {
        # s2 * corr_gauss_matrixvecC(x, y, theta)
        apply(x, 1, function(xx) {self$kone(xx, y, theta=theta, s2=s2)})
      } else if (is.matrix(y)) {
        # s2 * corr_gauss_matrixvecC(y, x, theta)
        apply(y, 1, function(yy) {self$kone(yy, x, theta=theta, s2=s2)})
      } else {
        self$kone(x, y, theta=theta, s2=s2)
      }
    },
    #' @description Calculate covariance between two points
    #' @param x vector.
    #' @param y vector.
    #' @param beta Correlation parameters. Log of theta.
    #' @param theta Correlation parameters.
    #' @param s2 Variance parameter.
    kone = function(x, y, beta, theta, s2) {
      # Kernels that inherit should implement this or k.
    },
    #' @description Starting point for parameters for optimization
    #' @param jitter Should there be a jitter?
    #' @param y Output
    #' @param beta_est Is beta being estimated?
    #' @param s2_est Is s2 being estimated?
    param_optim_start = function(jitter=F, y, beta_est=self$beta_est,
                                 s2_est=self$s2_est) {
      # Use current values for theta, partial MLE for s2
      # vec <- c(log(self$theta, 10),
      #  log(sum((y - mu) * solve(R, y - mu)) / n), 10)

      vec <- numeric(0)
      if (beta_est) {
        b <- c(self$beta)
        if (jitter) {
          b <- b + rnorm(length(b), 0, 1)
        }
        b <- pmin(pmax(b, self$beta_lower), self$beta_upper)
        vec <- c(vec, b)
      }
      if (s2_est) {
        b <- c(self$logs2)
        if (jitter) {
          b <- b + rnorm(length(b), 0, 1)
        }
        b <- pmin(pmax(b, self$logs2_lower), self$logs2_upper)
        vec <- c(vec, b)
      }
      vec
    },
    #' @description Starting point for parameters for optimization
    #' @param jitter Should there be a jitter?
    #' @param y Output
    #' @param beta_est Is beta being estimated?
    #' @param s2_est Is s2 being estimated?
    param_optim_start0 = function(jitter=F, y, beta_est=self$beta_est,
                                  s2_est=self$s2_est) {
      vec <- numeric(0)
      if (beta_est) {
        b <- rep(0, length(self$beta))
        if (jitter) {
          b <- b + rnorm(length(b), 0, 1)
        }
        b <- pmin(pmax(b, self$beta_lower), self$beta_upper)
        vec <- c(vec, b)
      }
      if (s2_est) {
        b <- 0
        if (jitter) {
          b <- b + rnorm(length(b), 0, 1)
        }
        b <- pmin(pmax(b, self$logs2_lower), self$logs2_upper)
        vec <- c(vec, b)
      }
      vec
    },
    #' @description Upper bounds of parameters for optimization
    #' @param p_est Is p being estimated?
    #' @param beta_est Is beta being estimated?
    #' @param s2_est Is s2 being estimated?
    param_optim_lower = function(beta_est=self$beta_est, s2_est=self$s2_est) {
      # c(self$beta_lower, self$logs2_lower)
      if (beta_est && s2_est) {
        c(self$beta_lower, self$logs2_lower)
      } else if (beta_est) {
        self$beta_lower
      } else if (s2_est) {
        self$logs2_lower
      } else {
        c()
      }
    },
    #' @description Upper bounds of parameters for optimization
    #' @param p_est Is p being estimated?
    #' @param beta_est Is beta being estimated?
    #' @param s2_est Is s2 being estimated?
    param_optim_upper = function(beta_est=self$beta_est, s2_est=self$s2_est) {
      # c(self$beta_upper, self$logs2_upper)
      if (beta_est && s2_est) {
        c(self$beta_upper, self$logs2_upper)
      } else if (beta_est) {
        self$beta_upper
      } else if (s2_est) {
        self$logs2_upper
      } else {
        c()
      }
    },
    #' @description Set parameters from optimization output
    #' @param optim_out Output from optimization
    #' @param beta_est Is beta being estimated?
    #' @param s2_est Is s2 being estimated?
    set_params_from_optim = function(optim_out, beta_est=self$beta_est,
                                     s2_est=self$s2_est) {
      loo <- length(optim_out)
      if (beta_est) {
        self$beta <- optim_out[1:(self$beta_length)]
      }
      if (s2_est) {
        self$logs2 <- optim_out[loo]
        self$s2 <- 10 ^ self$logs2
        if (self$logs2 >= self$logs2_upper) {
          message(paste0("* s2 is at maximum value after optimizing. ",
                         "Check the fit to see it this caused a bad fit. ",
                         "Consider changing the kernel s2_upper ",
                         " or rescaling the data."))
        }
        if (self$logs2 <= self$logs2_lower) {
          message(paste0("* s2 is at minimum value after optimizing. ",
                         "Check the fit to see it this caused a bad fit. ",
                         "Consider changing the kernel s2_lower ",
                         " or rescaling the data."))
        }
      }
    },
    # dC_dparams = function(params=NULL, C, X, C_nonug) {
    #   Kernels that inherit from this must implement this.
    # },
    #' @description Calculate covariance matrix and its derivative
    #'  with respect to parameters
    #' @param params Kernel parameters
    #' @param X matrix of points in rows
    #' @param nug Value of nugget
    C_dC_dparams = function(params=NULL, X, nug) {
      s2 <- self$s2_from_params(params)
      C_nonug <- self$k(x=X, params=params)
      C <- C_nonug + diag(s2*nug, nrow(X))
      dC_dparams <- self$dC_dparams(params=params, X=X, C_nonug=C_nonug,
                                    C=C, nug=nug)
      list(C=C, dC_dparams=dC_dparams)
    },
    #' @description Get s2 from params vector
    #' @param params parameter vector
    #' @param s2_est Is s2 being estimated?
    s2_from_params = function(params, s2_est=self$s2_est) {
      # 10 ^ params[length(params)]
      if (s2_est && !is.null(params)) { # Is last if in params
        10 ^ params[length(params)]
      } else { # Else it is just using set value, not being estimated
        self$s2
      }
    }
  ),
  private = list(

  )
)
