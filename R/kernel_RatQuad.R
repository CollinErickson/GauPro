#' Rational Quadratic Kernel R6 class
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
#' @field alpha alpha value (the exponent). Between 0 and 2.
#' @field logalpha Log of alpha
#' @field logalpha_lower Lower bound for log of alpha
#' @field logalpha_upper Upper bound for log of alpha
#' @field alpha_est Should alpha be estimated?
#' @examples
#' k1 <- RatQuad$new(beta=0, alpha=0)
RatQuad <- R6::R6Class(
  classname = "GauPro_kernel_RatQuad",
  inherit = GauPro_kernel_beta,
  public = list(
    # beta = NULL,
    # beta_lower = NULL,
    # beta_upper = NULL,
    # beta_length = NULL,
    # s2 = NULL, # variance coefficient to scale correlation matrix to covariance
    # logs2 = NULL,
    # logs2_lower = NULL,
    # logs2_upper = NULL,
    alpha = NULL,
    logalpha = NULL,
    logalpha_lower = NULL,
    logalpha_upper = NULL,
    alpha_est = NULL,
    #' @description Initialize kernel object
    #' @param beta Initial beta value
    #' @param alpha Initial alpha value
    #' @param s2 Initial variance
    #' @param D Number of input dimensions of data
    #' @param beta_lower Lower bound for beta
    #' @param beta_upper Upper bound for beta
    #' @param beta_est Should beta be estimated?
    #' @param alpha_lower Lower bound for alpha
    #' @param alpha_upper Upper bound for alpha
    #' @param alpha_est Should alpha be estimated?
    #' @param s2_lower Lower bound for s2
    #' @param s2_upper Upper bound for s2
    #' @param s2_est Should s2 be estimated?
    initialize = function(beta, alpha=1, s2=1, D,
                          beta_lower=-8, beta_upper=6, beta_est=TRUE,
                          alpha_lower=0, alpha_upper=Inf, alpha_est=TRUE,
                          s2_lower=1e-8, s2_upper=1e8, s2_est=TRUE
    ) {
      super$initialize(beta=beta, s2=s2, D=D, beta_lower=beta_lower,
                       beta_upper=beta_upper, beta_est=beta_est,
                       s2_lower=s2_lower,s2_upper=s2_upper, s2_est=s2_est)
      self$alpha <- alpha
      self$logalpha <- log(alpha, 10)
      self$logalpha_lower <- log(alpha_lower, 10)
      self$logalpha_upper <- log(alpha_upper, 10)
      self$alpha_est <- alpha_est

    },
    #' @description Calculate covariance between two points
    #' @param x vector.
    #' @param y vector, optional. If excluded, find correlation
    #' of x with itself.
    #' @param beta Correlation parameters.
    #' @param logalpha A correlation parameter
    #' @param s2 Variance parameter.
    #' @param params parameters to use instead of beta and s2.
    k = function(x, y=NULL, beta=self$beta, logalpha=self$logalpha, s2=self$s2, params=NULL) {#browser()
      if (!is.null(params)) {
        lenparams <- length(params)
        # beta <- params[1:(lenpar-2)]
        # logalpha <- params[lenpar-1]
        # logs2 <- params[lenpar]

        if (self$beta_est) {
          beta <- params[1:self$beta_length]
        } else {
          beta <- self$beta
        }
        if (self$alpha_est) {
          logalpha <- params[1 + as.integer(self$beta_est) * self$beta_length]
        } else {
          logalpha <- self$logalpha
        }
        if (self$s2_est) {
          logs2 <- params[lenparams]
        } else {
          logs2 <- self$logs2
        }

        s2 <- 10^logs2
      } else {#browser()
        if (is.null(beta)) {beta <- self$beta}
        if (is.null(logalpha)) {logalpha <- self$logalpha}
        if (is.null(s2)) {s2 <- self$s2}
      }
      theta <- 10^beta
      alpha <- 10^logalpha
      if (is.null(y)) {
        if (is.matrix(x)) {#browser()
          # cgmtry <- try(val <- s2 * corr_gauss_matrix_symC(x, theta))
          val <- outer(1:nrow(x), 1:nrow(x),
                       Vectorize(function(i,j){
                         self$kone(x[i,],x[j,],theta=theta, alpha=alpha, s2=s2)
                       }))
          # if (inherits(cgmtry,"try-error")) {browser()}
          return(val)
        } else {
          return(s2 * 1)
        }
      }
      if (is.matrix(x) & is.matrix(y)) {
        # s2 * corr_gauss_matrixC(x, y, theta)
        outer(1:nrow(x), 1:nrow(y), Vectorize(function(i,j){self$kone(x[i,],y[j,],theta=theta, alpha=alpha, s2=s2)}))
      } else if (is.matrix(x) & !is.matrix(y)) {
        # s2 * corr_gauss_matrixvecC(x, y, theta)
        apply(x, 1, function(xx) {self$kone(xx, y, theta=theta, alpha=alpha, s2=s2)})
      } else if (is.matrix(y)) {
        # s2 * corr_gauss_matrixvecC(y, x, theta)
        apply(y, 1, function(yy) {self$kone(yy, x, theta=theta, alpha=alpha, s2=s2)})
      } else {
        self$kone(x, y, theta=theta, alpha=alpha, s2=s2)
      }
    },
    #' @description Find covariance of two points
    #' @param x vector
    #' @param y vector
    #' @param beta correlation parameters on log scale
    #' @param alpha A correlation parameter
    #' @param theta correlation parameters on regular scale
    #' @param s2 Variance parameter
    kone = function(x, y, beta, theta, alpha, s2) {
      if (missing(theta)) {theta <- 10^beta}
      # t1 <- self$sqrt
      r2 <- sum(theta * (x-y)^2)
      s2 * (1 + r2 / alpha) ^ -alpha
    },
    #' @description Derivative of covariance with respect to parameters
    #' @param params Kernel parameters
    #' @param X matrix of points in rows
    #' @param C_nonug Covariance without nugget added to diagonal
    #' @param C Covariance with nugget
    #' @param nug Value of nugget
    dC_dparams = function(params=NULL, X, C_nonug, C, nug) {#browser(text = "Make sure all in one list")
      n <- nrow(X)
      # if (is.null(params)) {params <- c(self$beta, self$logalpha, self$logs2)}

      lenparams <- length(params)
      if (lenparams > 0) {
        if (self$beta_est) {
          beta <- params[1:self$beta_length]
        } else {
          beta <- self$beta
        }
        if (self$alpha_est) {
          logalpha <- params[1 + as.integer(self$beta_est) * self$beta_length]
        } else {
          logalpha <- self$logalpha
        }
        if (self$s2_est) {
          logs2 <- params[lenparams]
        } else {
          logs2 <- self$logs2
        }
      } else {
        beta <- self$beta
        logalpha <- self$logalpha
        logs2 <- self$logs2
      }

      # beta <- params[1:(lenparams - 2)]
      theta <- 10^beta
      # logalpha <- params[lenparams-1]
      alpha <- 10^logalpha
      log10 <- log(10)
      # logs2 <- params[lenparams]
      s2 <- 10 ^ logs2

      if (missing(C_nonug)) { # Assume C missing too, must have nug
        C_nonug <- self$k(x=X, params=params)
        C <- C_nonug + diag(nug*s2, nrow(C_nonug))
      }

      lenparams_D <- self$beta_length*self$beta_est + 1*self$alpha_est +self$s2_est
      dC_dparams <- array(dim=c(lenparams_D, n, n), data=0)
      if (self$s2_est) {
        dC_dparams[lenparams_D,,] <- C * log10
      }
      if (self$beta_est) {
        for (k in 1:length(beta)) {
          for (i in seq(1, n-1, 1)) {
            for (j in seq(i+1, n, 1)) {
              r2 <- sum(theta * (X[i,]-X[j,])^2)
              t1 <- 1 + r2 / alpha
              dC_dparams[k,i,j] <- -C_nonug[i,j] * (X[i,k] - X[j,k])^2  / t1 * theta[k] * log10   #s2 * (1+t1) * exp(-t1) *-dt1dbk + s2 * dt1dbk * exp(-t1)
              dC_dparams[k,j,i] <- dC_dparams[k,i,j]
            }
          }
          for (i in seq(1, n, 1)) { # Get diagonal set to zero
            dC_dparams[k,i,i] <- 0
          }
        }
      }
      if (self$alpha_est) {
        # Grad for alpha
        alpha_ind <- lenparams_D - self$s2_est
        for (i in seq(1, n-1, 1)) {
          for (j in seq(i+1, n, 1)) {
            r2 <- sum(theta * (X[i,]-X[j,])^2)
            t1 <- 1 + r2 / alpha
            dC_dparams[alpha_ind, i,j] <- C_nonug[i,j] * (- log(t1) + r2 / alpha / t1) * alpha * log10
            dC_dparams[alpha_ind, j,i] <- dC_dparams[alpha_ind, i,j]
          }
        }
        for (i in seq(1, n, 1)) {
          # dC_dparams[lenparams-1, i,i] <- 0
          dC_dparams[alpha_ind, i,i] <- 0
        }
      }
      return(dC_dparams)
    },
    #' @description Derivative of covariance with respect to X
    #' @param XX matrix of points
    #' @param X matrix of points to take derivative with respect to
    #' @param theta Correlation parameters
    #' @param beta log of theta
    #' @param alpha parameter
    #' @param s2 Variance parameter
    dC_dx = function(XX, X, theta, beta=self$beta, alpha=self$alpha, s2=self$s2) {#browser()
      if (missing(theta)) {theta <- 10^beta}
      # p <- 10 ^ logp
      # alpha <- 10 ^ logalpha
      if (!is.matrix(XX)) {stop()}
      d <- ncol(XX)
      if (ncol(X) != d) {stop()}
      n <- nrow(X)
      nn <- nrow(XX)
      dC_dx <- array(NA, dim=c(nn, d, n))
      for (i in 1:nn) {
        for (j in 1:d) {
          for (k in 1:n) {
            # r <- sqrt(sum(theta * (XX[i,] - X[k,]) ^ 2))
            r2 <- sum(theta * (XX[i,] - X[k,])^2)
            CC <- s2 * (1 + r2 / alpha) ^ -alpha
            dC_dx[i, j, k] <- CC * (-1) / (1 + r2 / alpha) * 2 * theta[j] * (XX[i, j]-X[k, j]) #) * p[j] #* (XX[i, j] - X[k, j])
          }
        }
      }
      dC_dx
    },
    #' @description Starting point for parameters for optimization
    #' @param jitter Should there be a jitter?
    #' @param y Output
    #' @param beta_est Is beta being estimated?
    #' @param alpha_est Is alpha being estimated?
    #' @param s2_est Is s2 being estimated?
    param_optim_start = function(jitter=F, y, beta_est=self$beta_est,
                                 alpha_est=self$alpha_est, s2_est=self$s2_est) {
      # Use current values for theta, partial MLE for s2
      # vec <- c(log(self$theta, 10), log(sum((y - mu) * solve(R, y - mu)) / n), 10)
      if (beta_est) {vec <- c(self$beta)} else {vec <- c()}
      if (alpha_est) {vec <- c(vec, self$logalpha)} else {}
      if (s2_est) {vec <- c(vec, self$logs2)} else {}
      if (jitter && beta_est) {
        # vec <- vec + c(self$beta_optim_jitter,  0)
        vec[1:length(self$beta)] = vec[1:length(self$beta)] + rnorm(length(self$beta), 0, 1)
      }
      vec
    },
    #' @description Starting point for parameters for optimization
    #' @param jitter Should there be a jitter?
    #' @param y Output
    #' @param beta_est Is beta being estimated?
    #' @param alpha_est Is alpha being estimated?
    #' @param s2_est Is s2 being estimated?
    param_optim_start0 = function(jitter=F, y, beta_est=self$beta_est,
                                  alpha_est=self$alpha_est, s2_est=self$s2_est) {
      # Use 0 for theta, partial MLE for s2
      # vec <- c(rep(0, length(self$theta)), log(sum((y - mu) * solve(R, y - mu)) / n), 10)
      if (beta_est) {vec <- rep(0, self$beta_length)} else {vec <- c()}
      if (alpha_est) {vec <- c(vec, 1)} else {}
      if (s2_est) {vec <- c(vec, 0)} else {}
      if (jitter && beta_est) {
        vec[1:length(self$beta)] = vec[1:length(self$beta)] + rnorm(length(self$beta), 0, 1)
      }
      vec
    },
    #' @description Lower bounds of parameters for optimization
    #' @param beta_est Is beta being estimated?
    #' @param alpha_est Is alpha being estimated?
    #' @param s2_est Is s2 being estimated?
    param_optim_lower = function(beta_est=self$beta_est,
                                 alpha_est=self$alpha_est,
                                 s2_est=self$s2_est) {
      # c(self$beta_lower, self$logs2_lower)
      if (beta_est) {vec <- c(self$beta_lower)} else {vec <- c()}
      if (alpha_est) {vec <- c(vec, self$logalpha_lower)} else {}
      if (s2_est) {vec <- c(vec, self$logs2_lower)} else {}
      vec
    },
    #' @description Upper bounds of parameters for optimization
    #' @param beta_est Is beta being estimated?
    #' @param alpha_est Is alpha being estimated?
    #' @param s2_est Is s2 being estimated?
    param_optim_upper = function(beta_est=self$beta_est,
                                 alpha_est=self$alpha_est,
                                 s2_est=self$s2_est) {
      # c(self$beta_upper, self$logs2_upper)
      if (beta_est) {vec <- c(self$beta_upper)} else {vec <- c()}
      if (alpha_est) {vec <- c(vec, self$logalpha_upper)} else {}
      if (s2_est) {vec <- c(vec, self$logs2_upper)} else {}
      vec
    },
    #' @description Set parameters from optimization output
    #' @param optim_out Output from optimization
    #' @param beta_est Is beta being estimated?
    #' @param alpha_est Is alpha being estimated?
    #' @param s2_est Is s2 being estimated?
    set_params_from_optim = function(optim_out, beta_est=self$beta_est,
                                     alpha_est=self$alpha_est, s2_est=self$s2_est) {
      loo <- length(optim_out)
      if (beta_est) {
        self$beta <- optim_out[1:(self$beta_length)]
      }
      if (alpha_est) {
        self$logalpha <- optim_out[(1 + beta_est * self$beta_length)]
        self$alpha <- 10 ^ self$logalpha
      }
      if (s2_est) {
        self$logs2 <- optim_out[loo]
        self$s2 <- 10 ^ self$logs2
      }
    }
  )
)
