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
    initialize = function(beta, s2=1, beta_lower=-8, beta_upper=6,
                          s2_lower=1e-8, s2_upper=1e8) {
      self$beta <- beta
      self$beta_length <- length(beta)
      # if (length(theta) == 1) {
      #   self$theta <- rep(theta, self$d)
      # }
      self$beta_lower <- beta_lower
      self$beta_upper <- beta_upper

      self$s2 <- s2
      self$logs2 <- log(s2, 10)
      self$logs2_lower <- log(s2_lower, 10)
      self$logs2_upper <- log(s2_upper, 10)
    },
    k = function(x, y=NULL, beta=self$beta, s2=self$s2, params=NULL) {#browser()
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
          # cgmtry <- try(val <- s2 * corr_gauss_matrix_symC(x, theta))
          val <- outer(1:nrow(x), 1:nrow(x), Vectorize(function(i,j){self$kone(x[i,],x[j,],theta=theta, s2=s2)}))
          # if (inherits(cgmtry,"try-error")) {browser()}
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
      s2 * exp(-sqrt(sum(theta * (x-y)^2)))
    },
    dC_dparams = function(params=NULL, C, X, C_nonug) {#browser(text = "Make sure all in one list")
      if (is.null(params)) {params <- c(self$beta, self$logs2)}
      lenparams <- length(params)
      beta <- params[1:(lenparams - 1)]
      theta <- 10^beta
      log10 <- log(10)
      logs2 <- params[lenparams]
      s2 <- 10 ^ logs2
      dC_dlogs2 <- C * log10 #/ s2 * s2 *
      dC_dbetas <- rep(list(C_nonug), length(beta))
      n <- nrow(X)
      for (k in 1:length(beta)) {
        for (i in seq(1, n-1, 1)) {
          for (j in seq(i+1, n, 1)) {
            dC_dbetas[[k]][i,j] <- -1 * dC_dbetas[[k]][i,j] * (X[i,k] - X[j,k])^2 * theta[k] * log10 * .5 / (-log(C[i,j]/s2))
            dC_dbetas[[k]][j,i] <- dC_dbetas[[k]][i,j]
          }
        }
        for (i in seq(1, n, 1)) { # Get diagonal set to zero
          dC_dbetas[[k]][i,i] <- 0
        }
        # Trying this, delete
        # dC_dbetas[[k]] <- -1 * dC_dbetas[[k]]
      }

      mats <- c(dC_dbetas, list(dC_dlogs2))
      return(list(dC_dparams=mats,
                  s2
      ))
    }
  )
)
