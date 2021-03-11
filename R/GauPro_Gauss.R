# correlation function should implement:
# corr: name of correlation
# corr_func
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



#' Corr Gauss GP using inherited optim
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @useDynLib GauPro
#' @importFrom Rcpp evalCpp
#' @importFrom stats optim
# @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' n <- 12
#' x <- matrix(seq(0,1,length.out = n), ncol=1)
#' y <- sin(2*pi*x) + rnorm(n,0,1e-1)
#' gp <- GauPro_Gauss$new(X=x, Z=y, parallel=FALSE)
GauPro_Gauss <- R6::R6Class(classname = "GauPro_Gauss",
     inherit = GauPro_base,
     public = list(
       corr = "Gauss",
       #corr_func = NULL,
       theta = NULL, # theta of length D, will repeat values if nonseparable
       theta_length = NULL, # number of unique theta's
       theta_map = NULL, # Vector of length D, ith component is which theta group ith dimension belongs to, eg c(1,2,1) means first and third dimensions share a theta
       theta_short = NULL, # vector of unique theta's, length D if separable, length 1 if nonsep. Use this for optimization, then transform to full theta using map
       #theta_short_length = NULL,
       separable = NULL,
       initialize = function(X, Z, verbose=0, separable=T, useC=F,useGrad=T,
                             parallel=FALSE,
                             nug=1e-6, nug.min=1e-8, nug.est=T,
                             param.est=T,
                             theta = NULL, theta_short = NULL, theta_map = NULL,
                             ...) {
         super$initialize(X=X,Z=Z,verbose=verbose,useC=useC,useGrad=useGrad,
                          parallel=parallel,
                          nug=nug, nug.min=nug.min, nug.est=nug.est, param.est=param.est)


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

         if (is.null(theta) & is.null(theta_short)) {
           self$theta_short <- rep(1, self$theta_length)
           self$theta <- self$theta_short[self$theta_map]
         } else if (!is.null(theta) & is.null(theta_short)) {
           if (!is.null(theta_map)) {
             stop("Don't give in theta and theta_map")
           }
           self$theta_short <- theta
           self$theta <- theta
         } else if (is.null(theta) & !is.null(theta_short)) {
           if (is.null(theta_map)) {
             stop("Don't give theta_short without giving theta_map")
           }
           self$theta_short <- theta_short
           self$theta <- self$theta_short[self$theta_map]
         } else if (!is.null(theta) & !is.null(theta_short)) {
           stop("Don't give both theta and theta_short")
         }



         self$fit()
         invisible(self)
       },
       #corr_func = GauPro::corr_gauss_matrix, # This is in R/corr, but copied here, DOESN'T WORK
       corr_func = function(x, x2=NULL, theta=self$theta) {
         if (is.null(x2)) corr_gauss_matrix_symC(x, theta)
         else corr_gauss_matrixC(x, x2, theta)
       },
       deviance_theta = function (theta) {
         self$deviance(theta=theta, nug=self$nug)
       },
       deviance_theta_log = function (beta) {
         theta <- 10^beta
         self$deviance_theta(theta=theta)
       },
       deviance = function(theta=self$theta, nug=self$nug) {
         if (length(theta) < self$D) {theta = theta[self$theta_map]} # if not fully separable, map out to full theta
         #Gaussian_devianceC(theta, nug, self$X, self$Z) # adding a try in case C can't chol()
         try.dev.C <- try(Gaussian_devianceC(theta, nug, self$X, self$Z))
         if (inherits(try.dev.C, "try-error")) {return(Inf)}
         try.dev.C
       },
       #deviance_out = function(...){Gaussian_devianceC(...)},
       deviance_grad = function(theta=NULL, nug=self$nug, joint=NULL, overwhat=if (self$nug.est) "joint" else "theta") {
         if (length(theta) < self$D) {theta = theta[self$theta_map]} # if not fully separable, map out to full theta
         gr <- Gaussian_deviance_gradC(theta=theta, nug=nug, X=self$X, Z=self$Z, overwhat=overwhat)
         if (!self$separable & overwhat!="nug") { # Map down if theta not full dimensional
           tgr <- sapply(1:self$theta_length, function(ii){sum(gr[which(self$theta_map == ii)])})
           if (overwhat == "joint") { # If joint, add nugget to end and return
             return(c(tgr, tail(gr,1)))
           }
           return(tgr) # else just return theta grad
         }
         gr
       },
       #deviance_grad_out = function(...){Gaussian_deviance_gradC(...)},
       deviance_fngr = function (theta=NULL, nug=NULL, overwhat=if (self$nug.est) "joint" else "theta") {#browser()
         if (length(theta) < self$D) {theta = theta[self$theta_map]} # if not fully separable, map out to full theta
         fngr <- Gaussian_deviance_fngrC(theta=theta, nug=nug, X=self$X, Z=self$Z, overwhat=overwhat)
         if (!self$separable & overwhat!="nug") {#browser() # Map down if theta not full dimensional
           gr <- fngr$gr
           tgr <- sapply(1:self$theta_length, function(ii){sum(gr[which(self$theta_map == ii)])})
           if (overwhat == "joint") { # If joint, add nugget to end and return
             return(list(fn=fngr$fn, gr=c(tgr, tail(gr,1))))
           }
           return(list(fn=fngr$fn, gr=tgr)) # else just return theta grad
         }
         fngr
       },
       #deviance_fngr_out = function(...){Gaussian_deviance_fngrC(...)},
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
       #deviance_grad = function(theta=NULL, nug=self$nug, joint=NULL, overwhat=if (self$nug.est) "joint" else "theta") {
       #  self$deviance_grad_out(theta=theta, nug=nug, X=self$X, Z=self$Z, overwhat=overwhat)
       #},
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
       #deviance_fngr = function (theta=NULL, nug=NULL, overwhat=if (self$nug.est) "joint" else "theta") {
       #  self$deviance_fngr_out(theta=theta, nug=nug, X=self$X, Z=self$Z, overwhat=overwhat)
       #},
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
         if (self$theta_length < self$D) {
           theta <- theta[self$theta_map]
           beta <- log(theta,10)
           joint <- c(beta, lognug)
         }
         tmp <- self$deviance_fngr(theta=theta, nug=nug, overwhat=overwhat)

         # Need to scale gradient, chain rule to get derivatives on log scale
         #if (length(tmp[[2]]) != length(joint)) {browser()}
         if (overwhat == "joint") {
           tmp[[2]] <- tmp[[2]] * 10^joint * log(10) # scale gradient only
         } else if (overwhat == "theta") {
           tmp[[2]] <- tmp[[2]] * theta * log(10) # scale gradient only
         } else if (overwhat == "nug") {
           tmp[[2]] <- tmp[[2]] * nug * log(10) # scale gradient only
         } else {
           print(paste("Overwhat is", overwhat, ", which is not accepted #523592"))
         }
         # Old below
         # tmp[[2]] <- tmp[[2]] * 10^joint * log(10) # scale gradient only

         tmp
       },


       get_optim_functions = function(param_update, nug.update){
         if (param_update & nug.update) {
           optim.func <- function(xx) {self$deviance_log2(joint=xx)}
           optim.grad <- function(xx) {self$deviance_log2_grad(joint=xx)}
           optim.fngr <- function(xx) {self$deviance_log2_fngr(joint=xx)}
         } else if (param_update & !nug.update) {
           optim.func <- function(xx) {self$deviance_log2(beta=xx)}
           optim.grad <- function(xx) {self$deviance_log2_grad(beta=xx, overwhat="theta")}
           optim.fngr <- function(xx) {self$deviance_log2_fngr(beta=xx, overwhat="theta")}
         } else if (!param_update & nug.update) {
           optim.func <- function(xx) {self$deviance_log2(lognug=xx)}
           optim.grad <- function(xx) {self$deviance_log2_grad(lognug=xx)}
           optim.fngr <- function(xx) {self$deviance_log2_fngr(lognug=xx)}
         } else {
           stop("Can't optimize over no variables")
         }
         return(list(optim.func=optim.func, optim.grad=optim.grad, optim.fngr=optim.fngr))
       },
       #optim.func = function(param_update, nug.update) {self$deviance},
       # optim.grad
       # optim.fngr - optimization fngr (list with function value and gradient vector including nugget)
       param_optim_lower = function() {# - lower bound of params
         rep(-5, self$theta_length)
       },
       param_optim_upper = function() {# - upper bound of params
         rep(5, self$theta_length)
       },
       param_optim_start = function() {# - current param values on optim scale
         pmin(pmax(log(self$theta_short, 10), -5), 5)
       },
       param_optim_start0 = function () {# - some central param values that can be used for optimization restarts
         rep(0, self$theta_length)
       },
       param_optim_jitter = function (param_value) { # only returns the jitter, not the param+jitter
         rnorm(self$theta_length,0,2)
       },

       update_params = function(restarts, param_update, nug.update) {
         if (param_update==FALSE && nug.update==FALSE) {
           return()
         }
         pars <- self$optim(
           restarts = restarts, param_update = param_update, nug.update = nug.update
         )$par

         # Check if nugget is below nug.min since lbfgs doesn't use bounds
         if (nug.update) {
           if (pars[length(pars)] < self$nug.min) {
             message("Below nug.min, setting nug to nug.min and reoptimizing #82387")
             self$nug <- self$nug.min
             nug.update=FALSE
             if (param_update) {
               # Set to best thetas and run from there
               self$theta_short <- 10 ^ pars[1:self$theta_length]
               self$theta <- self$theta_short[self$theta_map]

               # And rerun opt once fixing nugget
               pars <- self$optim(
                 restarts = 0, param_update = T, nug.update = F
               )$par
             }

           }
         }

         if (nug.update) {self$nug <- pars[length(pars)]}
         if (param_update) {
           self$theta_short <- 10 ^ pars[1:self$theta_length]
           self$theta <- self$theta_short[self$theta_map]
         }
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
       grad_dist = function (XX) {
         if (!is.matrix(XX)) {
           if (self$D == 1) XX <- matrix(XX, ncol=1)
           else if (length(XX) == self$D) XX <- matrix(XX, nrow=1)
           else stop('Predict input should be matrix')
         } else {
           if (ncol(XX) != self$D) {stop("Wrong dimension input")}
         }
         kx.xx <- self$corr_func(self$X, XX, theta=self$theta)

         # browser()
         xx <- as.numeric(XX)
         dkxx.x_dxx <- vapply(1:nrow(XX),
                              Vectorize(
                                function(k) {
                                  t(-2 * outer(1:self$N,
                                               1:self$D,
                                               Vectorize(
                                                 function(i,j) {
                                                   self$theta[j] * (XX[k, j] - self$X[i, j]) * kx.xx[i, k]}))
                                  )  #%*%self$Kinv %*% (self$Z - self$mu_hat)
                                }
                              )
                              , matrix(0, nrow=self$D, ncol=nrow(self$X))
         )
         cov_mat <- lapply(1:nrow(XX),
                function(i) {
                  (diag(2*self$theta) -
                    dkxx.x_dxx[,,i] %*% self$Kinv %*% t(dkxx.x_dxx[,,i])
                  ) * self$s2_hat
                  }
                )
         mean_vec <- lapply(1:nrow(XX),
                                function(i) {
                                  0 + dkxx.x_dxx[,,i] %*% self$Kinv %*% (self$Z - self$mu_hat)
                                }
                         )
         list(mean=mean_vec, cov=cov_mat)
       },
      hessian = function(XX, useC=self$useC) {#browser()

        if (!is.matrix(XX)) {
          if (self$D == 1) XX <- matrix(XX, ncol=1)
          else if (length(XX) == self$D) XX <- matrix(XX, nrow=1)
          else stop('Predict input should be matrix')
        } else {
          if (ncol(XX) != self$D) {stop("Wrong dimension input")}
        }
        #browser()
        hessian_func <- if (useC) {Gaussian_hessianC} else {Gaussian_hessianR}
        if (nrow(XX) == 1) {
          hess1 <- hessian_func(XX=XX, X=self$X, Z=self$Z, Kinv=self$Kinv, self$mu_hat, self$theta)
        } else { # more than one row
          hess1 <- lapply(1:nrow(XX), function(ii) {hessian_func(XX=XX[ii, ], X=self$X, Z=self$Z, Kinv=self$Kinv, self$mu_hat, self$theta)})
        }
        hess1
      },



       print = function() {
         cat("GauPro object of GauPro_Gauss\n")
         cat(paste0("\tD = ", self$D, ", N = ", self$N,"\n"))
         cat(paste0(c("\tTheta = ", signif(self$theta, 3), "\n")))
         cat(paste0("\tNugget = ", signif(self$nug, 3), "\n"))
         cat("\tRun update to add data and/or optimize again\n")
         cat("\tUse pred to get predictions at new points\n")
         invisible(self)
       }
     ),
     private = list(

     )
)
