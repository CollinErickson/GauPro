#' Exponential Kernel R6 class
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
#' k1 <- Exponential$new(beta=0)
Exponential <- R6::R6Class(classname = "GauPro_kernel_Exponential",
  inherit = GauPro_kernel_beta,
  public = list(
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
      if (is.null(y)) {
        if (is.matrix(x)) {
          val <- outer(1:nrow(x), 1:nrow(x), Vectorize(function(i,j){self$kone(x[i,],x[j,],theta=theta, s2=s2)}))
          return(val)
        } else {
          return(s2 * 1)
        }
      }
      if (is.matrix(x) & is.matrix(y)) {
        outer(1:nrow(x), 1:nrow(y), Vectorize(function(i,j){self$kone(x[i,],y[j,],theta=theta, s2=s2)}))
      } else if (is.matrix(x) & !is.matrix(y)) {
        apply(x, 1, function(xx) {self$kone(xx, y, theta=theta, s2=s2)})
      } else if (is.matrix(y)) {
        apply(y, 1, function(yy) {self$kone(yy, x, theta=theta, s2=s2)})
      } else {
        self$kone(x, y, theta=theta, s2=s2)
      }
    },
    kone = function(x, y, beta, theta, s2) {
      if (missing(theta)) {theta <- 10^beta}
      s2 * exp(-sqrt(sum(theta * (x-y)^2)))
    },
    dC_dparams = function(params=NULL, X, C_nonug, C, nug) {
      n <- nrow(X)
      if (is.null(params)) {params <- c(self$beta, self$logs2)}
      if (missing(C_nonug)) { # Assume C missing too, must have nug
        C_nonug <- self$k(x=X, params=params)
        C <- C_nonug + diag(nug*10^params[length(params)], nrow(C_nonug))
      }
      lenparams <- length(params)
      beta <- params[1:(lenparams - 1)]
      theta <- 10^beta
      log10 <- log(10)
      logs2 <- params[lenparams]
      s2 <- 10 ^ logs2
      dC_dparams <- array(dim=c(lenparams, n, n))
      dC_dparams[lenparams,,] <- C * log10 # Deriv for logs2

      # Derivs for beta
      for (i in seq(1, n-1, 1)) {
        for (j in seq(i+1, n, 1)) {
          t1 <- -1 * C_nonug[i,j] * log10 * .5 / (-log(C[i,j]/s2))
          for (k in 1:length(beta)) {
            dC_dparams[k,i,j] <- t1 * (X[i,k] - X[j,k])^2 * theta[k]
            dC_dparams[k,j,i] <- dC_dparams[k,i,j]
          }
        }
      }
      for (i in seq(1, n, 1)) { # Get diagonal set to zero
        for (k in 1:length(beta)) {
          dC_dparams[k,i,i] <- 0
        }
      }

      return(dC_dparams)
    },
    dC_dx = function(XX, X, theta, beta=self$beta, s2=self$s2) {#browser()
      if (missing(theta)) {theta <- 10^beta}
      if (!is.matrix(XX)) {stop()}
      d <- ncol(XX)
      if (ncol(X) != d) {stop()}
      n <- nrow(X)
      nn <- nrow(XX)
      dC_dx <- array(NA, dim=c(nn, d, n))
      for (i in 1:nn) {
        for (j in 1:d) {
          for (k in 1:n) {
            r <- sqrt(sum(theta * (XX[i,] - X[k,]) ^ 2))
            dC_dx[i, j, k] <- - theta[j] * (XX[i, j] - X[k, j]) * s2 * exp(-r) / r
          }
        }
      }
      dC_dx
    }
  )
)
