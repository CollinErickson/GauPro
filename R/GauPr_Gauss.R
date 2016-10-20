# correlation function should implement:
# corr: name of correlation
# corr_func
# deviance_out
# deviance
# deviance_grad
# deviance_fngr
# optim
# optimRestart
# optimBayes
# update_params
# grad

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
        initialize = function(X, Z, verbose=0, separable=T, useC=F,useGrad=T,
                              parallel=T, nug.est=T,
                              theta_map = NULL,
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
          Gaussian_devianceC(theta, nug, self$X, self$Z)
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
          #browser()
          tmp[[2]] <- tmp[[2]] * 10^joint * log(10) # scale gradient only
          tmp
        },


        optim = function (restarts = 5, theta.update = T, nug.update = self$nug.est, parallel=self$parallel, parallel.cores=self$parallel.cores) {
          # Does parallel
          # Joint MLE search with L-BFGS-B, with restarts
          if (theta.update & nug.update) {
            optim.func <- function(xx) {self$deviance_log2(joint=xx)}
            grad.func <- function(xx) {self$deviance_log2_grad(joint=xx)}
            optim.fngr <- function(xx) {self$deviance_log2_fngr(joint=xx)}
          } else if (theta.update & !nug.update) {
            optim.func <- function(xx) {self$deviance_log2(beta=xx)}
            grad.func <- function(xx) {self$deviance_log2_grad(beta=xx)}
            optim.fngr <- function(xx) {self$deviance_log2_fngr(beta=xx)}
          } else if (!theta.update & nug.update) {
            optim.func <- function(xx) {self$deviance_log2(lognug=xx)}
            grad.func <- function(xx) {self$deviance_log2_grad(lognug=xx)}
            optim.fngr <- function(xx) {self$deviance_log2_fngr(lognug=xx)}
          } else {
            stop("Can't optimize over no variables")
          }

          # This will make sure it at least can start
          # Run before it sets initial parameters
          try.devlog <- try(devlog <- self$deviance_log2(), silent = T)
          if (inherits(try.devlog, "try-error")) {
            warning("Current nugget doesn't work, increasing it #31973")
            self$update_K_and_estimates() # This will increase the nugget
            devlog <- self$deviance_log2()
          }

          lower <- c()
          upper <- c()
          start.par <- c()
          start.par0 <- c() # Some default params
          if (theta.update) {
            lower <- c(lower, rep(-5, self$theta_length))
            upper <- c(upper, rep(7, self$theta_length))
            start.par <- c(start.par, log(self$theta_short, 10))
            start.par0 <- c(start.par0, rep(0, self$theta_length))
          }
          if (nug.update) {
            lower <- c(lower, log(self$nug.min,10))
            upper <- c(upper, Inf)
            start.par <- c(start.par, log(self$nug,10))
            start.par0 <- c(start.par0, -6)
          }

          # Find best params with optimization, start with current params in case all give error
          # Current params
          best <- list(par=c(log(self$theta_short, 10), log(self$nug,10)), value = devlog)
          if (self$verbose >= 2) {cat("Optimizing\n");cat("\tInitial values:\n");print(best)}
          details <- data.frame(start=paste(c(self$theta_short,self$nug),collapse=","),end=NA,value=best$value,func_evals=1,grad_evals=NA,convergence=NA, message=NA, stringsAsFactors=F)


          # runs them in parallel, first starts from current, rest are jittered or random
          sys_name <- Sys.info()["sysname"]
          if (sys_name == "Windows" | !self$parallel) {
            # Trying this so it works on Windows
            restarts.out <- lapply( 1:(1+restarts), function(i){self$optimRestart2(start.par=start.par, start.par0=start.par0, theta.update=theta.update, nug.update=nug.update, optim.func=optim.func, grad.func=grad.func, optim.fngr=optim.fngr, lower=lower, upper=upper, jit=(i!=1))})#, mc.cores = parallel.cores)
          } else { # Mac/Unix
            restarts.out <- parallel::mclapply(1:(1+restarts), function(i){self$optimRestart2(start.par=start.par, start.par0=start.par0, theta.update=theta.update, nug.update=nug.update, optim.func=optim.func, grad.func=grad.func, optim.fngr=optim.fngr,lower=lower, upper=upper, jit=(i!=1))}, mc.cores = parallel.cores)
          }
          new.details <- t(sapply(restarts.out,function(dd){dd$deta}))
          vals <- sapply(restarts.out,
                         function(ii){
                           if (inherits(ii$current,"try-error")){Inf}
                           else ii$details$val
                         }
          )
          bestparallel <- which.min(vals) #which.min(new.details$value)
          if (restarts.out[[bestparallel]]$current$val < best$val) {
            best <- restarts.out[[bestparallel]]$current
          }
          details <- rbind(details, new.details)

          if (self$verbose >= 2) {print(details)}
          if (nug.update) best$par[length(best$par)] <- 10 ^ (best$par[length(best$par)])
          best
        },
        optimRestart = function (start.par, start.par0, theta.update, nug.update, optim.func, grad.func, optim.fngr, lower, upper, jit=T) {
          # FOR lognug RIGHT NOW, seems to be at least as fast, up to 5x on big data, many fewer func_evals
          #    still want to check if it is better or not
          #browser()
          if (runif(1) < .33 & jit) { # restart near some spot to avoid getting stuck in bad spot
            start.par.i <- start.par0
            #print("start at zero par")
          } else { # jitter from current params
            start.par.i <- start.par
          }
          if (jit) {
            if (theta.update) {start.par.i[1:self$theta_length] <- start.par.i[1:self$theta_length] + rnorm(self$theta_length,0,2)} # jitter betas
            if (nug.update) {start.par.i[length(start.par.i)] <- start.par.i[length(start.par.i)] + min(4, rexp(1,1))} # jitter nugget
          }
          if (self$verbose >= 2) {cat("\tRestart (parallel): starts pars =",start.par.i,"\n")}
          current <- try(
            if (self$useGrad) {
              if (is.null(optim.fngr)) {
                lbfgs::lbfgs(optim.func, grad.func, start.par.i, invisible=1)
              } else {
                lbfgs_share(optim.fngr, start.par.i, invisible=1) # 1.7x speedup uses grad_share
              }
            } else {
              optim(start.par.i, optim.func, method="L-BFGS-B", lower=lower, upper=upper, hessian=F)
            }
          )#;browser()
          if (!inherits(current, "try-error")) {#browser()
            if (is.character(current[[1]])) { # means cholesky failed probably
              details.new <- data.frame(start=paste(signif(start.par.i,3),collapse=","),end="cholesky error probably",value=Inf,func_evals=NA,grad_evals=NA,convergence=NA, message=current[1], stringsAsFactors=F)
            } else {
              if (self$useGrad) {current$counts <- c(NA,NA);if(is.null(current$message))current$message=NA}
              details.new <- data.frame(start=paste(signif(start.par.i,3),collapse=","),end=paste(signif(current$par,3),collapse=","),value=current$value,func_evals=current$counts[1],grad_evals=current$counts[2],convergence=current$convergence, message=current$message, row.names = NULL, stringsAsFactors=F)
            }
          } else{
            details.new <- data.frame(start=paste(signif(start.par.i,3),collapse=","),end="try-error",value=Inf,func_evals=NA,grad_evals=NA,convergence=NA, message=current[1], stringsAsFactors=F)
          }
          list(current=current, details=details.new)
        },
        optimBayes = function(theta.update = T, nug.update = F & self$nug.est) {
          lower <- c()
          upper <- c()
          start.par <- c()
          start.par0 <- c() # Some default params
          if (theta.update) {
            lower <- c(lower, rep(-5, self$theta_length))
            upper <- c(upper, rep(7, self$theta_length))
            start.par <- c(start.par, log(self$theta_short, 10))
            start.par0 <- c(start.par0, rep(0, self$theta_length))
          }
          if (nug.update) {
            lower <- c(lower, self$nug.min)
            upper <- c(upper, Inf)
            start.par <- c(start.par, self$nug)
            start.par0 <- c(start.par0, 1e-6)
          }
          # start with spread of points
          n0 <- 10
          d <- self$theta_length + as.numeric(nug.update)
          u <- matrix(start.par0,nrow=n0, ncol=d, byrow = T)
          for(i in 1:n0){u[i,] <- u[i,] + rnorm(d,0,2)}
          v <- apply(u,1, self$deviance_log)
          bgp <- GauPro$new(u,v)
          #curve(bgp$predict(x),-5,7);points(u,v)
          #curve(bgp$predict(x)+2*bgp$predict(x,se=T)$se,-5,7,add=T,col=2)
          #curve(bgp$predict(x)-2*bgp$predict(x,se=T)$se,-5,7,add=T,col=2)
          EI <- function(uu) {
            pr <- bgp$predict(uu,se=T)
            g <- (min(v) - pr$m) / pr$se
            pr$se * (g * pnorm(g) + dnorm(g))
          }
          for(i in 1:10) {
            #curve(EI, -5, 7);abline(v=u)
            if (F) {
              curve(bgp$predict(x),-5,7);points(u,v)
              curve(bgp$predict(x)+2*bgp$predict(x,se=T)$se,-5,7,add=T,col=2)
              curve(bgp$predict(x)-2*bgp$predict(x,se=T)$se,-5,7,add=T,col=2)
            }
            #apply(u, 1, function(uuu){optim(uuu, EI, lower=-5,upper=7)})
            #replicate(1e1, {optim(runif(1,-5,7), EI, lower=-5,upper=7,method="L-BFGS-B", control=list(fnscale=-1))$val})
            #tmp <- replicate(1e1, {optim(runif(1,-5,7), EI, lower=-5,upper=7,method="L-BFGS-B", control=list(fnscale=-1))[[c('par','val')]]})
            # parallel doesn't really speed up this line, already fast
            tmp <- replicate(10,{ta <- optim(runif(self$theta_length,-5,7), EI, lower=-5,upper=7,method="L-BFGS-B", control=list(fnscale=-1));c(ta$par,ta$val)})
            best <- tmp[,which.max(tmp[self$theta_length+1,])]
            bestu <- best[1:self$theta_length];#abline(v=bestu,col=2)
            bestv <- self$deviance_log(beta=bestu) # this might be slow, maybe do parallel and add a bunch of new points each time?
            bgp$update(Xnew = bestu, Znew = bestv, restarts = 0)
            u <- rbind(u, bestu)
            v <- c(v, bestv)
          }
          # optim from points already have, can pass grad #start at all midpoints?
          # take next point
          # start near that point?
          # repeat until limit reached
          bestu <- u[which.min(v),]
          optim(bestu, bgp$pred, lower=-5,upper=7,method="L-BFGS-B")$par
        },


        update_params = function(restarts, param_update, nug.update) {
          pars <- self$optim(
                    restarts = restarts, theta.update = param_update, nug.update = nug.update
                  )$par
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
