#' Matern 3/2 Kernel R6 class
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
#' @field sqrt3 Saved value of square root of 3
#' @examples
#' k1 <- Matern32$new(beta=0)
Matern32 <- R6::R6Class(classname = "GauPro_kernel_Matern32",
  inherit = GauPro_kernel_beta,
  public = list(
    sqrt3 = sqrt(3),
    #' @description Calculate covariance between two points
    #' @param x vector.
    #' @param y vector, optional. If excluded, find correlation
    #' of x with itself.
    #' @param beta Correlation parameters.
    #' @param s2 Variance parameter.
    #' @param params parameters to use instead of beta and s2.
    k = function(x, y=NULL, beta=self$beta, s2=self$s2, params=NULL) {
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
      } else {#browser()
        if (is.null(beta)) {beta <- self$beta}
        if (is.null(s2)) {s2 <- self$s2}
      }
      theta <- 10^beta
      if (is.null(y)) {
        if (is.matrix(x)) {
          # val <- outer(1:nrow(x), 1:nrow(x), Vectorize(function(i,j){self$kone(x[i,],x[j,],theta=theta, s2=s2)}))
          val <- s2 * corr_matern32_matrix_symC(x, theta)
          return(val)
        } else {
          return(s2 * 1)
        }
      }
      if (is.matrix(x) & is.matrix(y)) {
        # outer(1:nrow(x), 1:nrow(y), Vectorize(function(i,j){self$kone(x[i,],y[j,],theta=theta, s2=s2)}))
        s2 * corr_matern32_matrixC(x, y, theta)
      } else if (is.matrix(x) & !is.matrix(y)) {
        # apply(x, 1, function(xx) {self$kone(xx, y, theta=theta, s2=s2)})
        s2 * corr_matern32_matvecC(x, y, theta)
      } else if (is.matrix(y)) {
        # apply(y, 1, function(yy) {self$kone(yy, x, theta=theta, s2=s2)})
        s2 * corr_matern32_matvecC(y, x, theta)
      } else {
        self$kone(x, y, theta=theta, s2=s2)
      }
    },
    #' @description Find covariance of two points
    #' @param x vector
    #' @param y vector
    #' @param beta correlation parameters on log scale
    #' @param theta correlation parameters on regular scale
    #' @param s2 Variance parameter
    kone = function(x, y, beta, theta, s2) {
      if (missing(theta)) {theta <- 10^beta}
      r <- sqrt(sum(theta * (x-y)^2))
      t1 <- self$sqrt3 * r
      s2 * (1 + t1) * exp(-t1)
    },
    #' @description Derivative of covariance with respect to parameters
    #' @param params Kernel parameters
    #' @param X matrix of points in rows
    #' @param C_nonug Covariance without nugget added to diagonal
    #' @param C Covariance with nugget
    #' @param nug Value of nugget
    dC_dparams = function(params=NULL, X, C_nonug, C, nug) {
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

      lenparams_D <- self$beta_length*self$beta_est + self$s2_est
      dC_dparams <- array(dim=c(lenparams_D, n, n)) # Return as array
      if (self$s2_est) {
        dC_dparams[lenparams_D, , ] <- C * log10 # Deriv for logs2
      }

      # Loop to set beta derivatives
      if (self$beta_est) {
        for (i in seq(1, n-1, 1)) {
          for (j in seq(i+1, n, 1)) {
            tx2 <- sum(theta * (X[i,]-X[j,])^2)
            t1 <- sqrt(3 * tx2)
            t3 <- C[i,j] * (1/(1+t1) - 1) * self$sqrt3 * log10
            sqrttx2 <- sqrt(tx2)
            for (k in 1:length(beta)) {
              dt1dbk <- .5 * (X[i,k] - X[j,k])^2 / sqrttx2
              dC_dparams[k,i,j] <- t3 * dt1dbk * theta[k]   #s2 * (1+t1) * exp(-t1) *-dt1dbk + s2 * dt1dbk * exp(-t1)
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

      return(dC_dparams=dC_dparams)
    },
    #' @description Derivative of covariance with respect to X
    #' @param XX matrix of points
    #' @param X matrix of points to take derivative with respect to
    #' @param theta Correlation parameters
    #' @param beta log of theta
    #' @param s2 Variance parameter
    dC_dx = function(XX, X, theta, beta=self$beta, s2=self$s2) {
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
            dC_dx[i, j, k] <- -3 * s2 * r * exp(-self$sqrt3 * r) * theta[j] * (XX[i, j] - X[k, j]) / r
          }
        }
      }
      dC_dx
    }
  )
)
