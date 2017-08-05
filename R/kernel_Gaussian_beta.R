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



#' Gaussian Kernel R6 class
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
#' k1 <- Gaussian_beta$new(beta=0)
Gaussian_beta <- R6::R6Class(classname = "GauPro_kernel_Gaussian_beta",
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
    k = function(x, y=NULL, beta=self$beta, s2=self$s2, params=NULL) {
      if (!is.null(params)) {
        lenpar <- length(params)
        beta <- params[1:(lenpar-1)]
        logs2 <- params[lenpar]
        s2 <- 10^logs2
      } else {#browser()
        if (is.null(beta)) {beta <- self$beta}
        if (is.null(s2)) {s2 <- self$s2}
      }
      theta <- 10^beta
      if (is.null(y)) {
        if (is.matrix(x)) {#browser()
          cgmtry <- try(val <- s2 * corr_gauss_matrix_symC(x, theta))
          if (inherits(cgmtry,"try-error")) {browser()}
          return(val)
        } else {
          return(s2 * 1)
        }
      }
      if (is.matrix(x) & is.matrix(y)) {
        s2 * corr_gauss_matrixC(x, y, theta)
      } else if (is.matrix(x) & !is.matrix(y)) {
        s2 * corr_gauss_matrixvecC(x, y, theta)
      } else if (is.matrix(y)) {
        s2 * corr_gauss_matrixvecC(y, x, theta)
      } else {
        s2 * exp(-sum(theta * (x-y)^2))
      }
    },
    dC_dparams = function(params=NULL, X, C_nonug, C, nug) {#browser(text = "Make sure all in one list")
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
      dC_dparams[lenparams,,] <- C * log10 #/ s2 * s2 *
      # dC_dparams <- rep(list(C_nonug), length(beta))
      for (k in 1:length(beta)) {
        for (i in seq(1, n-1, 1)) {
          for (j in seq(i+1, n, 1)) {
            dC_dparams[k,i,j] <- - C_nonug[i,j] * (X[i,k] - X[j,k])^2 * theta[k] * log10
            dC_dparams[k,j,i] <- dC_dparams[k,i,j]
          }
        }
        for (i in seq(1, n, 1)) { # Get diagonal set to zero
          dC_dparams[k,i,i] <- 0
        }
      }

      # mats <- c(dC_dbetas, list(dC_dlogs2))
      return(dC_dparams)
    },
    C_dC_dparams = function(params=NULL, X, nug) {#browser(text = "Make sure all in one list")
      n <- nrow(X)
      if (is.null(params)) {params <- c(self$beta, self$logs2)}
      lenparams <- length(params)
      beta <- params[1:(lenparams - 1)]
      theta <- 10^beta
      log10 <- log(10)
      logs2 <- params[lenparams]
      s2 <- 10 ^ logs2

      # Calculate C
      C_nonug <- self$k(x=X, beta=beta, s2=s2)
      C <- C_nonug + diag(nug*s2, nrow(C_nonug))
      dC_dparams <- array(dim=c(lenparams, n, n))
      dC_dparams[lenparams,,] <- C * log10 #/ s2 * s2 *
      # dC_dbetas <- rep(list(C_nonug), length(beta))
      # n <- nrow(X)
      for (k in 1:length(beta)) {
        for (i in seq(1, n-1, 1)) {
          for (j in seq(i+1, n, 1)) {
            dC_dparams[k,i,j] <- - C[i,j] * (X[i,k] - X[j,k])^2 * theta[k] * log10
            dC_dparams[k,j,i] <- dC_dparams[k,i,j]
          }
        }
        for (i in seq(1, n, 1)) { # Get diagonal set to zero
          dC_dparams[k,i,i] <- 0
        }
      }

      # mats <- c(dC_dbetas, list(dC_dlogs2))
      return(list(C = C, dC_dparams))
    },
    dC_dx = functions(XX, X, theta, beta=self$beta, s2=self$s2) {
      if (missing(theta)) {theta <- 10^beta}
      if (!is.matrix(XX)) {stop()}
      d <- ncol(XX)
      if (ncol(X) != d) {stop()}
      n <- nrow(X)
      nn <- nrow(XX)
      dC_dx <- array(NA, dim=c(d, n, n))
      for (i in 1:nn) {
        for (j in 1:d) {
          for (k in 1:n) {
            dC_dx[i, j, k] <- -2 * theta[j] * (XX[i, j] - X[k, j]) * s2 * exp(-sum(theta * (XX[i,] - X[k,]) ^ 2))
          }
        }
      }
      dC_dx
    }
  )
)
