# Kernels should implement:
# k kernel function for two vectors
# update_params
# get_optim_functions: return optim.func, optim.grad, optim.fngr
# param_optim_lower - lower bound of params
# param_optim_upper - upper
# param_optim_start - current param values
# param_optim_start0 - some central param values that can be used for
#                      optimization restarts
# param_optim_jitter - how to jitter params in optimization

# Suggested
# deviance
# deviance_grad
# deviance_fngr
# grad



#' Gaussian Kernel R6 class
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
#' @examples
#' k1 <- Gaussian$new(beta=0)
Gaussian <- R6::R6Class(classname = "GauPro_kernel_Gaussian",
  inherit = GauPro_kernel_beta,
  public = list(
    # initialize = function(beta, s2=1, beta_lower=-8, beta_upper=6,
    #                       s2_lower=1e-8, s2_upper=1e8) {browser()
    #   self$beta <- beta
    #   self$beta_length <- length(beta)
    #   # if (length(theta) == 1) {
    #   #   self$theta <- rep(theta, self$d)
    #   # }
    #   self$beta_lower <- beta_lower
    #   self$beta_upper <- beta_upper
    #
    #   self$s2 <- s2
    #   self$logs2 <- log(s2, 10)
    #   self$logs2_lower <- log(s2_lower, 10)
    #   self$logs2_upper <- log(s2_upper, 10)
    # },
    #' @description Calculate covariance between two points
    #' @param x vector.
    #' @param y vector, optional. If excluded, find correlation
    #' of x with itself.
    #' @param beta Correlation parameters.
    #' @param s2 Variance parameter.
    #' @param params parameters to use instead of beta and s2.
    k = function(x, y=NULL, beta=self$beta, s2=self$s2, params=NULL) {
      if (!is.null(params)) {
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
        s2 <- 10 ^ logs2
      } else {#browser()
        if (is.null(beta)) {beta <- self$beta}
        if (is.null(s2)) {s2 <- self$s2}
      }
      theta <- 10^beta
      if (is.null(y)) {
        if (is.matrix(x)) {
          # cgmtry <- try(val <- s2 * corr_gauss_matrix_symC(x, theta))
          # if (inherits(cgmtry,"try-error")) {browser()}
          # return(val) # arma version isn't actually faster?
          return(s2 * corr_gauss_matrix_symC(x, theta))
          # return(s2 * corr_gauss_matrix_sym_armaC(x, theta))
        } else {
          return(s2 * 1)
        }
      }
      if (is.matrix(x) & is.matrix(y)) {
        s2 * corr_gauss_matrixC(x, y, theta)
        # Using Rcppparallel can be 2x faster but fails Travis
        # if (self$D >= 12 || nrow(x) < 30) {
        #   s2 * corr_gauss_matrixC(x, y, theta)
        # } else { # parallel only faster for small D and many rows
        #   s2 * corr_gauss_matrixCpar(x, y, theta)
        # }
        # s2 * corr_gauss_matrix_armaC(x, y, theta) # arma not actually faster?
        # corr_gauss_matrix_armaC(x, y, theta, s2)
      } else if (is.matrix(x) & !is.matrix(y)) {
        s2 * corr_gauss_matrixvecC(x, y, theta)
      } else if (is.matrix(y)) {
        s2 * corr_gauss_matrixvecC(y, x, theta)
      } else {
        s2 * exp(-sum(theta * (x-y)^2))
      }
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

      # if (inherits(try(diag(nug*s2, nrow(C_nonug))), "try-error")){browser()}
      # if (is.null(params)) {params <- c(self$beta, self$logs2)}
      if (missing(C_nonug)) { # Assume C missing too, must have nug
        C_nonug <- self$k(x=X, params=params)
        C <- C_nonug + diag(nug*s2, nrow(C_nonug))
      }

      lenparams_D <- self$beta_length*self$beta_est + self$s2_est

      # I wrote Rcpparmadillo function to speed this up a lot hopefully
      useR <- FALSE
      if (useR) {
        dC_dparams <- array(dim=c(lenparams_D, n, n), data=0)
        if (self$s2_est) {
          dC_dparams[lenparams_D,,] <- C * log10 #/ s2 * s2 *
        }
        # dC_dparams <- rep(list(C_nonug), length(beta))
        if (self$beta_est) {
          for (k in 1:length(beta)) {
            for (i in seq(1, n-1, 1)) {
              for (j in seq(i+1, n, 1)) {
                # if (inherits(try(C_nonug[i,j] * (X[i,k] - X[j,k])^2 *
                #           theta[k] * log10), "try-error")) {browser()}
                dC_dparams[k,i,j] <- - C_nonug[i,j] * (X[i,k] - X[j,k])^2 *
                                      theta[k] * log10
                dC_dparams[k,j,i] <- dC_dparams[k,i,j]
              }
            }
            for (i in seq(1, n, 1)) { # Get diagonal set to zero
              dC_dparams[k,i,i] <- 0
            }
          }
        }

      } else {
        dC_dparams <- kernel_gauss_dC(X, theta, C_nonug, self$s2_est,
                                      self$beta_est, lenparams_D, s2*nug)
      }
      # mats <- c(dC_dbetas, list(dC_dlogs2))
      return(dC_dparams)
    },
    #' @description Calculate covariance matrix and its derivative
    #'  with respect to parameters
    #' @param params Kernel parameters
    #' @param X matrix of points in rows
    #' @param nug Value of nugget
    C_dC_dparams = function(params=NULL, X, nug) {
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

      # if (is.null(params)) {params <- c(self$beta, self$logs2)}
      # beta <- params[1:(lenparams - 1)]
      theta <- 10^beta
      log10 <- log(10)
      # logs2 <- params[lenparams]
      s2 <- 10 ^ logs2
      # Calculate C
      C_nonug <- self$k(x=X, beta=beta, s2=s2)
      C <- C_nonug + diag(nug*s2, nrow(C_nonug))

      lenparams_D <- self$beta_length*self$beta_est + self$s2_est

      # I wrote Rcpparmadillo function to speed this up a lot hopefully
      useR <- FALSE
      if (useR) {
        dC_dparams <- array(dim=c(lenparams_D, n, n), data=0)
        if (self$s2_est) {
          dC_dparams[lenparams_D,,] <- C * log10 #/ s2 * s2 *
        }
        # dC_dbetas <- rep(list(C_nonug), length(beta))
        # n <- nrow(X)
        if (self$beta_est) {
          for (k in 1:length(beta)) {
            for (i in seq(1, n-1, 1)) {
              for (j in seq(i+1, n, 1)) {
                dC_dparams[k,i,j] <- - C[i,j] * (X[i,k] - X[j,k])^2 *
                                            theta[k] * log10
                dC_dparams[k,j,i] <- dC_dparams[k,i,j]
              }
            }
            for (i in seq(1, n, 1)) { # Get diagonal set to zero
              dC_dparams[k,i,i] <- 0
            }
          }
        }
      } else {
        dC_dparams <- kernel_gauss_dC(X, theta, C_nonug, self$s2_est,
                                      self$beta_est, lenparams_D, s2*nug)
      }
      # kernel_gauss_dC(X, theta, C_nonug, self$s2_est,
      #                  self$beta_est, lenparams_D, s2*nug)
      # mats <- c(dC_dbetas, list(dC_dlogs2))
      return(list(C = C, dC_dparams))
    },
    # dC_dx = function(XX, X, theta, beta=self$beta, s2=self$s2) {#browser()
    #   if (missing(theta)) {theta <- 10^beta}
    #   if (!is.matrix(XX)) {stop("XX must be matrix")}
    #   d <- ncol(XX)
    #   if (ncol(X) != d) {stop("XX and X must have same number")}
    #   n <- nrow(X)
    #   nn <- nrow(XX)
    #   dC_dx <- array(NA, dim=c(nn, d, n))
    #   for (i in 1:nn) {
    #     for (j in 1:d) {
    #       for (k in 1:n) {
    #         dC_dx[i, j, k] <- -2 * theta[j] * (XX[i, j] - X[k, j]) *
    #                             s2 * exp(-sum(theta * (XX[i,] - X[k,]) ^ 2))
    #       }
    #     }
    #   }
    #   dC_dx
    # },
    # Below is updated version using arma, was called dC_dx_arma before
    #' @description Derivative of covariance with respect to X
    #' @param XX matrix of points
    #' @param X matrix of points to take derivative with respect to
    #' @param theta Correlation parameters
    #' @param beta log of theta
    #' @param s2 Variance parameter
    dC_dx = function(XX, X, theta, beta=self$beta, s2=self$s2) {
      if (missing(theta)) {theta <- 10^beta}
      if (!is.matrix(XX)) {stop("XX must be matrix")}
      if (ncol(X) != ncol(XX)) {stop("XX and X must have same number")}
      corr_gauss_dCdX(XX, X, theta, s2)
    },
    #' @description Second derivative of covariance with respect to X
    #' @param XX matrix of points
    #' @param X matrix of points to take derivative with respect to
    #' @param theta Correlation parameters
    #' @param beta log of theta
    #' @param s2 Variance parameter
    d2C_dx2 = function(XX, X, theta, beta=self$beta, s2=self$s2) {
      if (missing(theta)) {theta <- 10^beta}
      if (!is.matrix(XX)) {stop("XX must be matrix")}
      d <- ncol(XX)
      if (ncol(X) != d) {stop("X and XX must have same # of columns")}
      n <- nrow(X)
      nn <- nrow(XX)
      d2C_dx2 <- array(NA, dim=c(nn, d, d, n))
      for (i in 1:nn) {
        for (k in 1:n) {
          Cik <- s2 * exp(-sum(theta * (XX[i,] - X[k,]) ^ 2))
          if (d > 1) {
            for (j1 in 1:(d-1)) {
              for (j2 in (j1+1):d) {
                d2C_dx2[i, j1, j2, k] <- 4 * theta[j1] *
                                        (XX[i, j1] - X[k, j1]) * theta[j2] *
                                        (XX[i, j2] - X[k, j2]) * Cik
                d2C_dx2[i, j2, j1, k] <- d2C_dx2[i, j1, j2, k]
              }
            }
          }
          for (j in 1:d) {
            d2C_dx2[i, j, j, k] <- -2 * theta[j] * Cik +
                                4 * theta[j]^2 * (XX[i, j] - X[k, j])^2 * Cik
          }
        }
      }
      d2C_dx2
    },
    #' @description Second derivative of covariance with respect to
    #' X and XX each once.
    #' @param XX matrix of points
    #' @param X matrix of points to take derivative with respect to
    #' @param theta Correlation parameters
    #' @param beta log of theta
    #' @param s2 Variance parameter
    d2C_dudv = function(XX, X, theta, beta=self$beta, s2=self$s2) {
      if (missing(theta)) {theta <- 10^beta}
      if (!is.matrix(XX)) {stop("XX must be matrix")}
      d <- ncol(XX)
      if (ncol(X) != d) {stop("X and XX must have same # of columns")}
      n <- nrow(X)
      nn <- nrow(XX)
      d2C_dx2 <- array(NA, dim=c(nn, d, d, n))
      for (i in 1:nn) {
        for (k in 1:n) {
          Cik <- s2 * exp(-sum(theta * (XX[i,] - X[k,]) ^ 2))
          if (d > 1) {
            for (j1 in 1:(d-1)) {
              for (j2 in (j1+1):d) {
                d2C_dx2[i, j1, j2, k] <- - 4 * theta[j1] *
                                        (XX[i, j1] - X[k, j1]) * theta[j2] *
                                        (XX[i, j2] - X[k, j2]) * Cik
                d2C_dx2[i, j2, j1, k] <- d2C_dx2[i, j1, j2, k]
              }
            }
          }
          for (j in 1:d) {
            d2C_dx2[i, j, j, k] <- 2 * theta[j] * Cik -
                                4 * theta[j]^2 * (XX[i, j] - X[k, j])^2 * Cik
          }
        }
      }
      d2C_dx2
    },
    #' @description Second derivative of covariance with respect to X and XX
    #' when they equal the same value
    #' @param XX matrix of points
    #' @param theta Correlation parameters
    #' @param beta log of theta
    #' @param s2 Variance parameter
    d2C_dudv_ueqvrows = function(XX, theta, beta=self$beta, s2=self$s2) {
      # Calculates derivative of C w.r.t. each component evaluated for
      #  both components equal to rows of XX
      # Vectorized version of d2C_dudv for u=v for rows of XX
      # Name is for "u equal v for rows of XX"
      # Much simpler since XX-X terms go to zero when XX=X
      # For m1 matrix, following two are equal, this version 2.5x faster
      # lapply(1:nrow(m1), function(i) {gp$kernel$d2C_dudv(XX = m1[i,,drop=F],
      #                                           X = m1[i,,drop=F])[1,,,1]})
      # gp$kernel$d2C_dudv_ueqvrows(XX = m1)
      if (missing(theta)) {theta <- 10^beta}
      if (!is.matrix(XX)) {stop("XX must be matrix")}
      d <- ncol(XX)
      nn <- nrow(XX)
      d2C_dx2 <- array(0, dim=c(nn, d, d))
      for (j in 1:d) {
        # Not multiplied by C since C=1 when u=v
        d2C_dx2[, j, j] <- 2 * theta[j] * s2
      }
      d2C_dx2
    }
  )
)
