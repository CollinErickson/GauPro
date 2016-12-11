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
#' @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' n <- 12
#' x <- matrix(seq(0,1,length.out = n), ncol=1)
#' y <- sin(2*pi*x) + rnorm(n,0,1e-1)
#' gp <- GauPro(X=x, Z=y, parallel=FALSE)
GauPr_LB <- R6::R6Class(classname = "GauPr_LB",
     inherit = GauPr,
     public = list(
       corr = "Lifted Brownian",
       #corr_func = NULL,
       theta = NULL, # theta of length D, will repeat values if nonseparable
       theta_length = NULL, # number of unique theta's
       theta_map = NULL, # Vector of length D, ith component is which theta group ith dimension belongs to, eg c(1,2,1) means first and third dimensions share a theta
       theta_short = NULL, # vector of unique theta's, length D if separable, length 1 if nonsep. Use this for optimization, then transform to full theta using map
       theta_min = 1e-5,
       theta_max = 1e5,
       beth = NULL, # in (0,1)
       beth_min = 1e-8,
       beth_max = 1-1e-8,
       gam = NULL, # in (0,1]
       gam_min = 1e-8,
       gam_max = 1,
       separable = NULL,
       initialize = function(X, Z, verbose=0, separable=T, useC=F,useGrad=T,
                             parallel=T, nug.est=T,
                             theta_map = NULL,
                             beth=0.5, gam=1,
                             ...) {
         super$initialize(X=X,Z=Z,verbose=verbose,useC=useC,useGrad=useGrad,
                          parallel=parallel, nug.est=nug.est)


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
         self$beth <- beth
         self$gam <- gam


         self$fit()
         invisible(self)
       },
       #corr_func = GauPro::corr_gauss_matrix, # This is in R/corr, but copied here, DOESN'T WORK
       corr_func = function(x, x2=NULL, theta=self$theta) {
         #if (is.null(x2)) corr_gauss_matrix_symC(x, theta)
         #else corr_gauss_matrixC(x, x2, theta)
       },
       psi = function(x, theta=self$theta, beth=self$beth, gam=self$gam) {
         (1 + sum(a * x ^ 2) ^ gam) ^ beth
       },
       cov_func = function(x, x2=NULL, theta=self$theta, beth=self$beth, gam=self$gam) {
         if   (is.null(x2)) {
           psi(x, theta=theta, beth=beth, gam=gam) +
             psi(x2, theta=theta, beth=beth, gam=gam) +
             psi(x-x2, theta=theta, beth=beth, gam=gam) - 1
         } else {
           2 * psi(x, theta=theta, beth=beth, gam=gam) +
             psi(x-x2, theta=theta, beth=beth, gam=gam) - 1
         }
       },
       deviance_theta = function (theta) {
         self$deviance(theta=theta, beth=self$beth, gam=self$gam, nug=self$nug)
       },
       deviance_theta_log = function (beta) {
         theta <- 10^beta
         self$deviance_theta(theta=theta)
       },
       deviance = function(theta=self$theta, beth=self$beth, gam=self$gam, nug=self$nug) {
         if (length(theta) < self$D) {theta = theta[self$theta_map]} # if not fully separable, map out to full theta
         #Gaussian_devianceC(theta, nug, self$X, self$Z)
         Sigma <- self$cov_func(x=self$X, theta=theta, beth=beth, gam=gam) + nug * diag(self$N)
         log(det(Sigma)) + t(self$Z) %*% solve(Sigma, self$Z)
       },
       #deviance_out = function(...){Gaussian_devianceC(...)},
       dpsi_dtheta = function(x,theta=self$theta,beth=self$beth,gam=self$gam) {
         magha2 <- sum(theta*x^2)^gam
         beth * (1 + magha2)^(beth-1) * gam * (magha2)^(gam-1) * x^2
       },
       dpsi_dbeth = function(x,theta=self$theta,beth=self$beth,gam=self$gam) {
         psi(x=x,theta=theta, beth=beth, gam=gam) * log(1+sum(theta*x^2)^gam)
       },
       dpsi_dgam = function(x,theta=self$theta,beth=self$beth,gam=self$gam) {
         #beth * psi(x=x,theta=theta, beth=beth - 1, gam=gam) * sum(theta*x^2)^gam * log(sum(theta*x^2)^gam)
         magha2 <- sum(theta*x^2)^gam # Protect against log(0)
         magha2gam <- magha2^gam
         if (magha2 == 0) {0}
         else {beth * (1+magha2gam)^(beth-1) * maghagam * log(magha2)}
       },
       dc_dtheta = function(x1, x2, theta=self$theta,beth=self$beth,gam=self$gam) {
         dpsi_dtheta(x=x1,theta=theta,beth=beth,gam=gam) +
           dpsi_dtheta(x=x2,theta=theta,beth=beth,gam=gam) +
           dpsi_dtheta(x=x1-x2,theta=theta,beth=beth,gam=gam)
       },
       dc_dbeth = function(x1, x2, theta=self$theta,beth=self$beth,gam=self$gam) {
         dpsi_dbeth(x=x1,theta=theta,beth=beth,gam=gam) +
           dpsi_dbeth(x=x2,theta=theta,beth=beth,gam=gam) +
           dpsi_dbeth(x=x1-x2,theta=theta,beth=beth,gam=gam)
       },
       dc_dgam = function(x1, x2, theta=self$theta,beth=self$beth,gam=self$gam) {
         dpsi_dgam(x=x1,theta=theta,beth=beth,gam=gam) +
           dpsi_dgam(x=x2,theta=theta,beth=beth,gam=gam) +
           dpsi_dgam(x=x1-x2,theta=theta,beth=beth,gam=gam)
       },
       deviance_gradR = function(theta=self$theta, beth=self$beth, gam=self$gam, nug=self$nug) {
         N <- self$N
         dSigma_dbeth = matrix(NA, N, N)
         dSigma_dgam = matrix(NA, N, N)
         dSigma_dtheta = replicate(self$D, matrix(NA, N, N))
         for(i in 1:N) {
           for(j in 1:N) {
             dSigma_dbeth <- dc_dbeth(self$X[i,], self$X[j,], theta=theta,beth=beth,gam=gam)
             dSigma_dgam <- dc_dgam(self$X[i,], self$X[j,], theta=theta,beth=beth,gam=gam)
             dSigma_dtheta[i,j,] <- dc_dtheta(self$X[i,], self$X[j,], theta=theta,beth=beth,gam=gam)
           }
         }
         Sigma = self$cov_func(theta=theta,beth=beth,gam=gam) + nug * diag(N)
         Sigmainv_y <- solve(Sigma, self$Z)
         dD_dbeth <- sum(diag(solve(Sigma, dSigma_dbeth))) + t(Sigmainv_y) %*% dSigma_dbeth %*% dSigma_dbeth
         dD_dgam <- sum(diag(solve(Sigma, dSigma_dgam))) + t(Sigmainv_y) %*% dSigma_dgam %*% dSigma_dgam
         dD_dtheta <- numeric(N)
         for(i in 1:N) {
           dD_dtheta[i] <- sum(diag(solve(Sigma, dSigma_dtheta[,,i]))) + t(Sigmainv_y) %*% dSigma_dtheta[,,i] %*% dSigma_dgam
         }
         dD_dnug <- sum(diag(solve(Sigma))) + sum(Sigmainv_y^2)
         c(dD_dtheta, dD_dbeth, dD_dtheta, dD_dnug)
       },
       deviance_grad = function(theta=NULL, nug=self$nug, joint=NULL, overwhat=if (self$nug.est) "joint" else "theta") {
         if (length(theta) < self$D) {theta = theta[self$theta_map]} # if not fully separable, map out to full theta
         gr <- self$deviance_gradR(theta=theta, beth=beth, gam=gam, nug=nug)#, X=self$X, Z=self$Z, overwhat=overwhat)
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
         #browser()
         tmp[[2]] <- tmp[[2]] * 10^joint * log(10) # scale gradient only
         tmp
       },


       get_optim_functions = function(param_update, nug.update){
         if (param_update & nug.update) {
           optim.func <- function(xx) {self$deviance_log2(joint=xx)}
           optim.grad <- function(xx) {self$deviance_log2_grad(joint=xx)}
           optim.fngr <- function(xx) {self$deviance_log2_fngr(joint=xx)}
         } else if (param_update & !nug.update) {
           optim.func <- function(xx) {self$deviance_log2(beta=xx)}
           optim.grad <- function(xx) {self$deviance_log2_grad(beta=xx)}
           optim.fngr <- function(xx) {self$deviance_log2_fngr(beta=xx)}
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
         c(rep(self$beta_min, self$theta_length),
           self$beth_min,
           self$gam_min
         )
       },
       param_optim_upper = function() {# - upper bound of params
         c(rep(self$beta_max, self$theta_length),
           self$beth_max,
           self$gam_max
         )
       },
       param_optim_start = function() {# - current param values on optim scale
         c(log(self$theta_short, 10),
           self$beth,
           self$gam
         )
       },
       param_optim_start0 = function () {# - some central param values that can be used for optimization restarts
         c(rep(0, self$theta_length),
           .5,
           .99
         )
       },
       param_optim_jittered = function (param_value) { # returns the param+jitter
         beth <- param_value[self$theta_length + 1]
         gam <- param_value[self$theta_length + 2]
         c(param_value[1:self$theta_length] + rnorm(self$theta_length,0,2),
           beth, # not jittering for now, too difficult
           gam
         )
       },

       update_params = function(restarts, param_update, nug.update) {
         pars <- self$optim(
           restarts = restarts, param_update = param_update, nug.update = nug.update
         )$par
         if (nug.update) {self$nug <- pars[length(pars)]}
         if (param_update) {
           self$theta_short <- 10 ^ pars[1:self$theta_length]
           self$theta <- self$theta_short[self$theta_map]
           self$beth <- pars[self$theta_length + 1]
           self$gam <- pars[self$theta_length + 2]
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




       print = function() {
         cat("GauPro object\n")
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
