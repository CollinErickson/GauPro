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
    initialize = function(p, alpha, s2=1,
                          p_lower=0, p_upper=1e2, p_est=TRUE,
                          alpha_lower=0, alpha_upper=1e2, alpha_est=TRUE,
                          s2_lower=1e-8, s2_upper=1e8, s2_est=TRUE
    ) {
      self$alpha <- alpha
      self$logalpha <- log(alpha, 10)
      self$logalpha_lower <- log(alpha_lower, 10)
      self$logalpha_upper <- log(alpha_upper, 10)
      self$alpha_est <- alpha_est
      self$p <- p
      self$p_length <- length(p)
      self$logp <- log(p, 10)
      self$logp_lower <- log(p_lower, 10)
      self$logp_upper <- log(p_upper, 10)
      self$p_est <- p_est
      self$s2 <- s2
      self$logs2 <- log(s2, 10)
      self$logs2_lower <- log(s2_lower, 10)
      self$logs2_upper <- log(s2_upper, 10)
      self$s2_est <- s2_est

    },
    k = function(x, y=NULL, logp=self$logp, logalpha=self$logalpha, s2=self$s2, params=NULL) {#browser()
      if (!is.null(params)) {
        lenpar <- length(params)
        logp <- params[1:(lenpar-2)]
        logalpha <- params[lenpar-1]
        logs2 <- params[lenpar]
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
    kone = function(x, y, logp, p, alpha, s2) {
      if (missing(p)) {p <- 10^logp}
      out <- s2 * exp(-sum(alpha*sin(p * (x-y))^2))
      if (any(is.nan(out))) {browser()}
      out
    },
    dC_dparams = function(params=NULL, X, C_nonug, C, nug) {#browser(text = "Make sure all in one list")
      n <- nrow(X)
      if (is.null(params)) {params <- c(self$logp, self$logalpha, self$logs2)}
      if (missing(C_nonug)) { # Assume C missing too, must have nug
        C_nonug <- self$k(x=X, params=params)
        C <- C_nonug + diag(nug*10^params[length(params)], nrow(C_nonug))
      }
      lenparams <- length(params)
      logp <- params[1:(lenparams - 2)]
      p <- 10^logp
      logalpha <- params[lenparams-1]
      alpha <- 10^logalpha
      log10 <- log(10)
      logs2 <- params[lenparams]
      s2 <- 10 ^ logs2
      dC_dparams <- array(dim=c(lenparams, n, n))
      dC_dparams[lenparams,,] <- C * log10
      for (k in 1:length(logp)) {
        for (i in seq(1, n-1, 1)) {
          for (j in seq(i+1, n, 1)) {
            r2 <- sum(p * (X[i,]-X[j,])^2)
            dC_dparams[k,i,j] <- -C[i,j] * alpha * sin(2*p[k]*(X[i,k] - X[j,k])) * (X[i,k] - X[j,k]) * p[k] * log10
            dC_dparams[k,j,i] <- dC_dparams[k,i,j]
          }
        }
        for (i in seq(1, n, 1)) { # Get diagonal set to zero
          dC_dparams[k,i,i] <- 0
        }
      }
      # Grad for logalpha
      for (i in seq(1, n-1, 1)) {
        for (j in seq(i+1, n, 1)) {
          r2 <- -sum(sin(p * (X[i,]-X[j,]))^2)
          dC_dparams[lenparams-1, i,j] <- C[i,j] * r2 * alpha * log10
          dC_dparams[lenparams-1, j,i] <- dC_dparams[lenparams-1, i,j]
        }
      }
      for (i in seq(1, n, 1)) {
        dC_dparams[lenparams-1, i,i] <- 0
      }
      return(dC_dparams)
    },
    C_dC_dparams = function(params=NULL, X, nug) {
      s2 <- self$s2_from_params(params)
      C_nonug <- self$k(x=X, params=params)
      C <- C_nonug + diag(s2*nug, nrow(X))
      dC_dparams <- self$dC_dparams(params=params, X=X, C_nonug=C_nonug, C=C, nug=nug)
      list(C=C, dC_dparams=dC_dparams)
    },
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
    param_optim_lower = function(p_est=self$p_est,
                                 alpha_est=self$alpha_est,
                                 s2_est=self$s2_est) {
      # c(self$logp_lower, self$logs2_lower)
      if (p_est) {vec <- c(self$logp_lower)} else {vec <- c()}
      if (alpha_est) {vec <- c(vec, self$logalpha_lower)} else {}
      if (s2_est) {vec <- c(vec, self$logs2_lower)} else {}
      vec
    },
    param_optim_upper = function(p_est=self$p_est,
                                 alpha_est=self$alpha_est,
                                 s2_est=self$s2_est) {
      # c(self$logp_upper, self$logs2_upper)
      if (p_est) {vec <- c(self$logp_upper)} else {vec <- c()}
      if (alpha_est) {vec <- c(vec, self$logalpha_upper)} else {}
      if (s2_est) {vec <- c(vec, self$logs2_upper)} else {}
      vec
    },
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
