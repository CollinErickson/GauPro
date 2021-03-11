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



#' Periodic Kernel R6 class
#'
#' p is the period for each dimension, a is a single number for scaling
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
#' @field p Parameter for correlation
#' @field p_est Should p be estimated?
#' @field logp Log of p
#' @field logp_lower Lower bound of logp
#' @field logp_upper Upper bound of logp
#' @field p_length length of p
#' @field alpha Parameter for correlation
#' @field alpha_est Should alpha be estimated?
#' @field logalpha Log of alpha
#' @field logalpha_lower Lower bound of logalpha
#' @field logalpha_upper Upper bound of logalpha

#' @field s2 variance
#' @field s2_est Is s2 estimated?
#' @field logs2 Log of s2
#' @field logs2_lower Lower bound of logs2
#' @field logs2_upper Upper bound of logs2
#' @examples
#' k1 <- Periodic$new(p=1, alpha=1)
Periodic <- R6::R6Class(
  classname = "GauPro_kernel_Periodic",
  inherit = GauPro_kernel,
  public = list(
    p = NULL, # Period
    p_est = NULL,
    logp = NULL,
    logp_lower = NULL,
    logp_upper = NULL,
    p_length = NULL,
    s2 = NULL, # variance coefficient to scale correlation matrix to covariance
    s2_est = NULL,
    logs2 = NULL,
    logs2_lower = NULL,
    logs2_upper = NULL,
    alpha = NULL,
    logalpha = NULL,
    logalpha_lower = NULL,
    logalpha_upper = NULL,
    alpha_est = NULL,
    #' @description Initialize kernel object
    #' @param p Periodic parameter
    #' @param alpha Periodic parameter
    #' @param s2 Initial variance
    #' @param D Number of input dimensions of data
    #' @param p_lower Lower bound for p
    #' @param p_upper Upper bound for p
    #' @param p_est Should p be estimated?
    #' @param alpha_lower Lower bound for alpha
    #' @param alpha_upper Upper bound for alpha
    #' @param alpha_est Should alpha be estimated?
    #' @param s2_lower Lower bound for s2
    #' @param s2_upper Upper bound for s2
    #' @param s2_est Should s2 be estimated?
    initialize = function(p, alpha=1, s2=1, D,
                          p_lower=0, p_upper=1e2, p_est=TRUE,
                          alpha_lower=0, alpha_upper=1e2, alpha_est=TRUE,
                          s2_lower=1e-8, s2_upper=1e8, s2_est=TRUE
    ) {

      # Check p and D
      missing_p <- missing(p)
      missing_D    <- missing(D)
      if (missing_p && missing_D) {stop("Must give kernel p or D")}
      else if (missing_p) {p <- rep(pi, D)}
      else if (missing_D) {D <- length(p)}
      else {if (length(p) != D) {stop("p and D should have same length")}}

      self$D <- D


      self$alpha <- alpha
      self$logalpha <- log(alpha, 10)
      self$logalpha_lower <- log(alpha_lower, 10)
      self$logalpha_upper <- log(alpha_upper, 10)
      self$alpha_est <- alpha_est
      self$p <- p
      self$p_length <- length(p)
      self$logp <- log(p, 10)

      # Now set upper and lower so they have correct length
      # self$logp_lower <- log(p_lower, 10)
      # self$logp_upper <- log(p_upper, 10)
      # Setting logp_lower so dimensions are right
      logp_lower <- log(p_lower, 10)
      logp_upper <- log(p_upper, 10)
      self$logp_lower <- if (length(logp_lower) == self$p_length) {logp_lower}
      else if (length(logp_lower)==1) {rep(logp_lower, self$p_length)}
      else {stop("Error for kernel_Periodic logp_lower")}
      self$logp_upper <- if (length(logp_upper) == self$p_length) {logp_upper}
      else if (length(logp_upper)==1) {rep(logp_upper, self$p_length)}
      else {stop("Error for kernel_Periodic logp_upper")}

      self$p_est <- p_est
      self$s2 <- s2
      self$logs2 <- log(s2, 10)
      self$logs2_lower <- log(s2_lower, 10)
      self$logs2_upper <- log(s2_upper, 10)
      self$s2_est <- s2_est

    },
    #' @description Calculate covariance between two points
    #' @param x vector.
    #' @param y vector, optional. If excluded, find correlation
    #' of x with itself.
    #' @param logp Correlation parameters.
    #' @param logalpha Correlation parameters.
    #' @param s2 Variance parameter.
    #' @param params parameters to use instead of beta and s2.
    k = function(x, y=NULL, logp=self$logp, logalpha=self$logalpha, s2=self$s2, params=NULL) {#browser()
      if (!is.null(params)) {
        lenparams <- length(params)
        # logp <- params[1:(lenpar-2)]
        # logalpha <- params[lenpar-1]
        # logs2 <- params[lenpar]

        if (self$p_est) {
          logp <- params[1:self$p_length]
        } else {
          logp <- self$logp
        }
        if (self$alpha_est) {
          logalpha <- params[1 + as.integer(self$p_est) * self$p_length]
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
        if (is.null(logp)) {logp <- self$logp}
        if (is.null(logalpha)) {logalpha <- self$logalpha}
        if (is.null(s2)) {s2 <- self$s2}
      }
      p <- 10^logp
      alpha <- 10^logalpha
      if (is.null(y)) {
        if (is.matrix(x)) {#browser()
          val <- outer(1:nrow(x), 1:nrow(x),
                       Vectorize(function(i,j){
                         self$kone(x[i,],x[j,],p=p, alpha=alpha, s2=s2)
                       }))
          # if (inherits(cgmtry,"try-error")) {browser()}
          return(val)
        } else {
          return(s2 * 1)
        }
      }
      if (is.matrix(x) & is.matrix(y)) {
        outer(1:nrow(x), 1:nrow(y), Vectorize(function(i,j){self$kone(x[i,],y[j,],p=p, alpha=alpha, s2=s2)}))
      } else if (is.matrix(x) & !is.matrix(y)) {
        apply(x, 1, function(xx) {self$kone(xx, y, p=p, alpha=alpha, s2=s2)})
      } else if (is.matrix(y)) {
        apply(y, 1, function(yy) {self$kone(yy, x, p=p, alpha=alpha, s2=s2)})
      } else {
        self$kone(x, y, p=p, alpha=alpha, s2=s2)
      }
    },
    #' @description Find covariance of two points
    #' @param x vector
    #' @param y vector
    #' @param logp correlation parameters on log scale
    #' @param p correlation parameters on regular scale
    #' @param alpha correlation parameter
    #' @param s2 Variance parameter
    kone = function(x, y, logp, p, alpha, s2) {
      if (missing(p)) {p <- 10^logp}
      out <- s2 * exp(-sum(alpha*sin(p * (x-y))^2))
      if (any(is.nan(out))) {browser()}
      out
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

      if (lenparams > 0) {
        if (self$p_est) {
          logp <- params[1:self$p_length]
        } else {
          logp <- self$logp
        }
        if (self$alpha_est) {
          logalpha <- params[1 + as.integer(self$p_est) * self$p_length]
        } else {
          logalpha <- self$logalpha
        }
        if (self$s2_est) {
          logs2 <- params[lenparams]
        } else {
          logs2 <- self$logs2
        }
      } else {
        logp <- self$logp
        logalpha <- self$logalpha
        logs2 <- self$logs2
      }

      # lenparams <- length(params)
      # logp <- params[1:(lenparams - 2)]
      p <- 10^logp
      # logalpha <- params[lenparams-1]
      alpha <- 10^logalpha
      log10 <- log(10)
      # logs2 <- params[lenparams]
      s2 <- 10 ^ logs2

      # if (is.null(params)) {params <- c(self$logp, self$logalpha, self$logs2)}
      if (missing(C_nonug)) { # Assume C missing too, must have nug
        C_nonug <- self$k(x=X, params=params)
        C <- C_nonug + diag(nug*s2, nrow(C_nonug))
      }

      lenparams_D <- self$p_length*self$p_est + 1*self$alpha_est +self$s2_est
      dC_dparams <- array(dim=c(lenparams_D, n, n), data=0)
      if (self$s2_est) {
        dC_dparams[lenparams_D,,] <- C * log10
      }
      if (self$p_est) {
        for (k in 1:length(logp)) {
          for (i in seq(1, n-1, 1)) {
            for (j in seq(i+1, n, 1)) {
              r2 <- sum(p * (X[i,]-X[j,])^2)
              dC_dparams[k,i,j] <- -C_nonug[i,j] * alpha * sin(2*p[k]*(X[i,k] - X[j,k])) * (X[i,k] - X[j,k]) * p[k] * log10
              dC_dparams[k,j,i] <- dC_dparams[k,i,j]
            }
          }
          for (i in seq(1, n, 1)) { # Get diagonal set to zero
            dC_dparams[k,i,i] <- 0
          }
        }
      }
      # Grad for logalpha
      if (self$alpha_est) {
        alph_ind <- lenparams_D - as.integer(self$s2_est)
        for (i in seq(1, n-1, 1)) {
          for (j in seq(i+1, n, 1)) {
            r2 <- -sum(sin(p * (X[i,]-X[j,]))^2)
            dC_dparams[alph_ind, i,j] <- C_nonug[i,j] * r2 * alpha * log10
            dC_dparams[alph_ind, j,i] <- dC_dparams[alph_ind, i,j]
          }
        }
        for (i in seq(1, n, 1)) {
          dC_dparams[alph_ind, i,i] <- 0
        }
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
    #' @param logp log of p
    #' @param logalpha log of alpha
    #' @param s2 Variance parameter
    dC_dx = function(XX, X, logp=self$logp, logalpha=self$logalpha, s2=self$s2) {#browser()
      # if (missing(theta)) {theta <- 10^beta}
      p <- 10 ^ logp
      alpha <- 10 ^ logalpha
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
            CC <- s2 * exp(-sum(alpha * sin(p * (XX[i, ]-X[k, ]))^2))
            dC_dx[i, j, k] <- CC * (-alpha) * sin(2*p[j]*(XX[i, j]-X[k, j])) * p[j] #* (XX[i, j] - X[k, j])
          }
        }
      }
      dC_dx
    },
    #' @description Starting point for parameters for optimization
    #' @param jitter Should there be a jitter?
    #' @param y Output
    #' @param p_est Is p being estimated?
    #' @param alpha_est Is alpha being estimated?
    #' @param s2_est Is s2 being estimated?
    param_optim_start = function(jitter=F, y, p_est=self$p_est,
                                 alpha_est=self$alpha_est, s2_est=self$s2_est) {
      # Use current values for theta, partial MLE for s2
      # vec <- c(log(self$theta, 10), log(sum((y - mu) * solve(R, y - mu)) / n), 10)
      if (p_est) {vec <- c(self$logp)} else {vec <- c()}
      if (alpha_est) {vec <- c(vec, self$logalpha)} else {}
      if (s2_est) {vec <- c(vec, self$logs2)} else {}
      if (jitter && p_est) {
        # vec <- vec + c(self$logp_optim_jitter,  0)
        vec[1:length(self$logp)] = vec[1:length(self$logp)] + rnorm(length(self$logp), 0, 1)
      }
      vec
    },
    #' @description Starting point for parameters for optimization
    #' @param jitter Should there be a jitter?
    #' @param y Output
    #' @param p_est Is p being estimated?
    #' @param alpha_est Is alpha being estimated?
    #' @param s2_est Is s2 being estimated?
    param_optim_start0 = function(jitter=F, y, p_est=self$p_est,
                                  alpha_est=self$alpha_est, s2_est=self$s2_est) {
      # Use 0 for theta, partial MLE for s2
      # vec <- c(rep(0, length(self$theta)), log(sum((y - mu) * solve(R, y - mu)) / n), 10)
      if (p_est) {vec <- rep(0, self$p_length)} else {vec <- c()}
      if (alpha_est) {vec <- c(vec, 1)} else {}
      if (s2_est) {vec <- c(vec, 0)} else {}
      if (jitter && p_est) {
        vec[1:length(self$logp)] = vec[1:length(self$logp)] + rnorm(length(self$logp), 0, 1)
      }
      vec
    },
    #' @description Lower bounds of parameters for optimization
    #' @param p_est Is p being estimated?
    #' @param alpha_est Is alpha being estimated?
    #' @param s2_est Is s2 being estimated?
    param_optim_lower = function(p_est=self$p_est,
                                 alpha_est=self$alpha_est,
                                 s2_est=self$s2_est) {
      # c(self$logp_lower, self$logs2_lower)
      if (p_est) {vec <- c(self$logp_lower)} else {vec <- c()}
      if (alpha_est) {vec <- c(vec, self$logalpha_lower)} else {}
      if (s2_est) {vec <- c(vec, self$logs2_lower)} else {}
      vec
    },
    #' @description Upper bounds of parameters for optimization
    #' @param p_est Is p being estimated?
    #' @param alpha_est Is alpha being estimated?
    #' @param s2_est Is s2 being estimated?
    param_optim_upper = function(p_est=self$p_est,
                                 alpha_est=self$alpha_est,
                                 s2_est=self$s2_est) {
      # c(self$logp_upper, self$logs2_upper)
      if (p_est) {vec <- c(self$logp_upper)} else {vec <- c()}
      if (alpha_est) {vec <- c(vec, self$logalpha_upper)} else {}
      if (s2_est) {vec <- c(vec, self$logs2_upper)} else {}
      vec
    },
    #' @description Set parameters from optimization output
    #' @param optim_out Output from optimization
    #' @param p_est Is p being estimated?
    #' @param alpha_est Is alpha being estimated?
    #' @param s2_est Is s2 being estimated?
    set_params_from_optim = function(optim_out, p_est=self$p_est,
                                     alpha_est=self$alpha_est, s2_est=self$s2_est) {
      loo <- length(optim_out)
      if (p_est) {
        self$logp <- optim_out[1:(self$p_length)]
        self$p <- 10 ^ self$logp
      }
      if (alpha_est) {
        self$logalpha <- optim_out[(1 + p_est * self$p_length)]
        self$alpha <- 10 ^ self$logalpha
      }
      if (s2_est) {
        self$logs2 <- optim_out[loo]
        self$s2 <- 10 ^ self$logs2
      }
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
  )
)
