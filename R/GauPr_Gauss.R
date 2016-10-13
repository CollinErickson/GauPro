# correlation function should implement:
# corr: name of correlation
# corr_func
# deviance_out

# Should only require??
# deviance
# deviance_grad
# deviance_fngr
# optim_objective ???? to replace deviance_log2_fngr


#' Corr Gauss GP
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @useDynLib GauPro
#' @importFrom Rcpp evalCpp
#' @importFrom stats optim
#' @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' n <- 12
#' x <- matrix(seq(0,1,length.out = n), ncol=1)
#' y <- sin(2*pi*x) + rnorm(n,0,1e-1)
#' gp <- GauPro$new(X=x, Z=y, parallel=FALSE)
GauPr_Gauss <- R6::R6Class(classname = "GauPr_Gauss",
      inherit = GauPr,
      public = list(
        corr = "Gauss",
        #corr_func = NULL,
        theta = NULL, # theta of length D, will repeat values if nonseparable
        theta_length = NULL, # number of unique theta's
        theta_map = NULL, # Vector of length D, ith component is which theta group ith dimension belongs to, eg c(1,2,1) means first and third dimensions share a theta
        theta_short = NULL, # vector of unique theta's, length D if separable, length 1 if nonsep. Use this for optimization, then transform to full theta using map
        #theta_short_length = NULL,
        separable = NULL,
        mu_hat = NULL,
        s2_hat = NULL,
        corr = "Gauss",
        initialize = function(X, Z, verbose=0, separable=T, useC=F,useGrad=T,
                              parallel=T, useOptim2=T, nug.est=T,
                              theta_map = NULL,
                              ...) {
          self$X <- X
          self$Z <- matrix(Z, ncol=1)
          self$verbose <- verbose
          if (!is.matrix(self$X)) {
            if (length(self$X) == length(self$Z)) {
              self$X <- matrix(X, ncol=1)
            } else {
              stop("X and Z don't match")
            }
          }
          self$N <- nrow(self$X)
          self$D <- ncol(self$X)
          self$separable <- separable
          if (!is.null(theta_map)) {
            self$theta_map <- theta_map
            self$theta_length <- length(unique(theta_map))
          } else if (self$separable) {
            self$theta_map <- 1:self$D
            self$theta_length <- self$D
          } else{
            self$theta_map <- rep(1, self$D)
            self$theta_length <- 1
          }
          self$theta_short <- rep(1, self$theta_length)
          self$theta <- self$theta_short[self$theta_map]

          self$nug.est <- nug.est
          self$useC <- useC
          self$useGrad <- useGrad
          self$parallel <- parallel
          if (self$parallel) {self$parallel.cores <- parallel::detectCores()}
          else {self$parallel.cores <- 1}

          self$useOptim2 <- useOptim2

          self$fit()
          invisible(self)
        },
        simple = function() {print(133)},
        #corr_func = GauPro::corr_gauss_matrix, # This is in R/corr, but copied here, DOESN'T WORK
        corr_func = function(x, x2=NULL, theta) {
          if (is.null(x2)) corr_gauss_matrix_symC(x, theta)
          else corr_gauss_matrixC(x, x2, theta)
        },
        deviance_out = function(...){Gaussian_devianceC(...)},
        deviance_grad_out = function(...){Gaussian_deviance_gradC(...)},
        deviance_fngr_out = function(...){Gaussian_deviance_fngrC(...)},
        deviance_log = function (beta=NULL, nug=self$nug, joint=NULL) {
          if (!is.null(joint)) {
            beta <- joint[-length(joint)]
            theta <- 10^beta
            nug <- joint[length(joint)]
          } else {
            if (is.null(beta)) theta <- self$theta
            else theta <- 10^beta
          }
          self$deviance(theta=theta, nug=nug)
        },
        deviance_log2 = function (beta=NULL, lognug=NULL, joint=NULL) {  # joint deviance
          # This takes nug on log scale
          if (!is.null(joint)) {
            beta <- joint[-length(joint)]
            theta <- 10^beta
            nug <- 10^joint[length(joint)]
          } else {
            if (is.null(beta)) theta <- self$theta
            else theta <- 10^beta
            if (is.null(lognug)) nug <- self$nug
            else nug <- 10^lognug
          }
          self$deviance(theta=theta, nug=nug)
        },
        deviance_grad = function(theta=NULL, nug=self$nug, joint=NULL, overwhat=if (self$nug.est) "joint" else "theta") {
          self$deviance_grad_out(theta=theta, nug=nug, X=self$X, Z=self$Z, overwhat=overwhat)
        },
        deviance_log_grad = function (beta=NULL, nug=self$nug, joint=NULL, overwhat=if (self$nug.est) "joint" else "theta") {
          if (!is.null(joint)) {
            beta <- joint[-length(joint)]
            theta <- 10^beta
            nug <- joint[length(joint)]
          } else {
            if (is.null(beta)) {theta <- self$theta; beta <- log(theta, 10)}
            else theta <- 10^beta
          }
          dg <- self$deviance_gradC(theta=theta, nug=nug, overwhat=overwhat)
          if (overwhat != "nug") {dg[1:self$theta_length] <- dg[1:self$theta_length] * 10^beta * log(10)}
          dg
        },
        deviance_log2_grad = function (beta=NULL, lognug=NULL, joint=NULL, overwhat=if (self$nug.est) "joint" else "theta") {
          if (!is.null(joint)) {
            beta <- joint[-length(joint)]
            theta <- 10^beta
            lognug <- joint[length(joint)]
            nug <- 10^lognug
          } else {
            if (is.null(beta)) {theta <- self$theta; beta <- log(theta, 10)}
            else theta <- 10^beta
            if (is.null(lognug)) {nug <- self$nug; lognug <- log(nug, 10)}
            else nug <- 10^lognug
            joint <- c(beta, lognug)
          }
          self$deviance_gradC(theta=theta, nug=nug, overwhat=overwhat) * 10^joint * log(10)
        },
        deviance_fngr = function (theta=NULL, nug=NULL, overwhat=if (self$nug.est) "joint" else "theta") {
          self$deviance_fngr_out(theta=theta, nug=nug, X=self$X, Z=self$Z, overwhat=overwhat)
        },
        deviance_log2_fngr = function (beta=NULL, lognug=NULL, joint=NULL, overwhat=if (self$nug.est) "joint" else "theta") {
          if (!is.null(joint)) {
            beta <- joint[-length(joint)]
            theta <- 10^beta
            lognug <- joint[length(joint)]
            nug <- 10^lognug
          } else {
            if (is.null(beta)) {theta <- self$theta; beta <- log(theta, 10)}
            else theta <- 10^beta
            if (is.null(lognug)) nug <- self$nug
            else nug <- 10^lognug
            joint <- c(beta, lognug)
          }
          tmp <- self$deviance_fngr(theta=theta, nug=nug, overwhat=overwhat)
          tmp[[2]] <- tmp[[2]] * 10^joint * log(10) # scale gradient only
          tmp
        },
        grad = function (XX) {
          if (!is.matrix(XX)) {
            if (self$D == 1) XX <- matrix(XX, ncol=1)
            else if (length(XX) == self$D) XX <- matrix(XX, nrow=1)
            else stop('Predict input should be matrix')
          } else {
            if (ncol(XX) != self$D) {stop("Wrong dimension input")}
          }
          kx.xx <- self$corr_func(self$X, XX, theta=self$theta)

          grad1 <-   vapply(1:nrow(XX),
                            Vectorize(
                              function(k) {
                                t(-2 * outer(1:self$N, 1:self$D, Vectorize(function(i,j) {self$theta[j] * (XX[k, j] - self$X[i, j]) * kx.xx[i, k]}))
                                )  %*%self$Kinv %*% (self$Z - self$mu_hat)
                              }
                            )
                            , numeric(self$D)
          )
          if (self$D == 1) return(grad1)
          t(grad1)
        },
        grad_norm = function (XX) {
          grad1 <- self$grad(XX)
          if (!is.matrix(grad1)) return(abs(grad1))
          apply(grad1,1, function(xx) {sqrt(sum(xx^2))})
        },
        print = function() {
          cat("GauPro object\n")
          cat(paste0("\tD = ", self$D, ", N = ", self$N,"\n"))
          cat(paste0(c("\tTheta = ", signif(self$theta, 3), "\n")))
          cat(paste0("\tNugget = ", signif(self$nug, 3), "\n"))
          cat("\tRun update to add data and/or optimize again\n")
          cat("\tUse pred to get predictions at new points\n")
          invisible(self)
        }
        ,plot = plot
      ),
      private = list(

      )
)
