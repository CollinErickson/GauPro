#' Matern 5/2 Kernel R6 class
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
#' k1 <- Matern52$new(beta=0)
Matern52 <- R6::R6Class(classname = "GauPro_kernel_Matern52",
  inherit = GauPro_kernel_beta,
  public = list(
    sqrt5 = sqrt(5),
    k = function(x, y=NULL, beta=self$beta, s2=self$s2, params=NULL) {#browser()
      if (!is.null(params)) {
        # lenpar <- length(params)
        # beta <- params[1:(lenpar-1)]
        # logs2 <- params[lenpar]

        lenparams <- length(params)
        if (self$beta_est) {
          beta <- params[1:self$beta_length]
        } else {
          beta <- self$beta
        }
        if (self$s2_est) {
          logs2 <- params[lenparams]
        } else {
          logs2 <- self$logs2
        }

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
        # s2 * corr_gauss_matrixC(x, y, theta)
        outer(1:nrow(x), 1:nrow(y), Vectorize(function(i,j){self$kone(x[i,],y[j,],theta=theta, s2=s2)}))
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
    kone = function(x, y, beta, theta, s2) {
      if (missing(theta)) {theta <- 10^beta}
      r <- sqrt(sum(theta * (x-y)^2))
      t1 <- self$sqrt5 * r
      s2 * (1 + t1 + t1^2 / 3) * exp(-t1)
    },
    dC_dparams = function(params=NULL, X, C_nonug, C, nug) {#browser(text = "Make sure all in one list")
      n <- nrow(X)

      lenparams <- length(params)
      if (lenparams > 0) {
        if (self$beta_est) {
          beta <- params[1:self$beta_length]
        } else {
          beta <- self$beta
        }
        if (self$s2_est) {
          logs2 <- params[lenparams]
        } else {
          logs2 <- self$logs2
        }
      } else {
        beta <- self$beta
        logs2 <- self$logs2
      }

      # lenparams <- length(params)
      # beta <- params[1:(lenparams - 1)]
      theta <- 10^beta
      log10 <- log(10)
      # logs2 <- params[lenparams]
      s2 <- 10 ^ logs2

      # if (is.null(params)) {params <- c(self$beta, self$logs2)}
      if (missing(C_nonug)) { # Assume C missing too, must have nug
        C_nonug <- self$k(x=X, params=params)
        C <- C_nonug + diag(nug*s2, nrow(C_nonug))
      }
      dC_dparams <- array(dim=c(lenparams, n, n), data = 0)
      if (self$s2_est) {
        dC_dparams[lenparams,,] <- C * log10 # Deriv for logs2
      }

      # Deriv for beta
      if (self$beta_est) {
        for (i in seq(1, n-1, 1)) {
          for (j in seq(i+1, n, 1)) {
            tx2 <- sum(theta * (X[i,]-X[j,])^2)
            t1 <- sqrt(5 * tx2)
            t3 <- C[i,j] * ((1+2*t1/3)/(1+t1+t1^2/3) - 1) * self$sqrt5 * log10
            half_over_sqrttx2 <- .5 / sqrt(tx2)
            for (k in 1:length(beta)) {
              dt1dbk <- half_over_sqrttx2 * (X[i,k] - X[j,k])^2
              dC_dparams[k,i,j] <- t3 * dt1dbk * theta[k]
              dC_dparams[k,j,i] <- dC_dparams[k,i,j]
            }
          }
        }
        for (i in seq(1, n, 1)) { # Get diagonal set to zero
          for (k in 1:length(beta)) {
            dC_dparams[k,i,i] <- 0
          }
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
            dC_dx[i, j, k] <- (-5*r/3 - 5/3*self$sqrt5*r^2) * s2 * exp(-self$sqrt5 * r) * theta[j] * (XX[i, j] - X[k, j]) / r
          }
        }
      }
      dC_dx
    }
  )
)
