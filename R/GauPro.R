library(R6)
GauPro <- R6Class(classname = "GauPro",
  public = list(
    X = NULL,
    Z = NULL,
    N = NULL,
    D = NULL,
    corr = "Gauss",
    corr_func = NULL,
    nug = 1e-6,
    nug.min = 1e-8,
    nug.est = T,
    theta = NULL,
    theta_length = NULL,
    separable = NULL,
    mu_hat = NULL,
    s2_hat = NULL,
    K = NULL,
    Kchol = NULL,
    Kinv = NULL,
    verbose = 0,
    useC = TRUE,
    parallel = FALSE,
    parallel.cores = NULL,
    useOptim2 = F,
    initialize = function(X, Z, corr="Gauss", verbose=0, separable=T, useC=F,
                          parallel=T, useOptim2=F, ...) {#browser()
      for (item in list(...)) {
        self$add(item)
      }
      self$X <- X
      self$Z <- matrix(Z, ncol=1)
      #if (!is.null(corr)) {self$corr <- corr}
      self$corr <- corr
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
      if (self$separable) {
        self$theta_length <- self$D
        self$theta <- rep(1, self$theta_length)
      } else {
        self$theta_length <- 1
        self$theta <- 1
      }
      if (self$corr == "Gauss") {
        self$corr_func <- if (self$useC) {GauPro::corr_gauss_matrix}
                          else {GauPro::corr_gauss_matrix_noC}
      } else {
        stop("corr not specified or recognized")
      }

      self$useC <- useC
      self$parallel <- parallel
      if (self$parallel) {self$parallel.cores <- parallel::detectCores()}
      else {self$parallel.cores <- 1}

      self$useOptim2 <- useOptim2

      self$fit()
      invisible(self)
    },
    fit = function(X, Z) {
      self$update()
    },
    update_params = function () {#browser()
      self$K <- self$corr_func(self$X, theta=self$theta) + diag(self$nug, self$N)
      self$Kchol <- chol(self$K)
      self$Kinv <- chol2inv(self$Kchol)
      self$mu_hat <- sum(self$Kinv %*% self$Z) / sum(self$Kinv)
      self$s2_hat <- c(t(self$Z - self$mu_hat) %*% self$Kinv %*% (self$Z - self$mu_hat) / self$N)
    },
    predict = function(XX, se.fit=F, covmat=F) {
      self$pred(XX=XX, se.fit=se.fit, covmat=covmat)
    },
    pred = function(XX, se.fit=F, covmat=F) {#browser()
      if (!is.matrix(XX)) {
        if (self$D == 1) XX <- matrix(XX, ncol=1)
        else if (length(XX) == self$D) XX <- matrix(XX, nrow=1)
        else stop('Predict input should be matrix')
      }
      #covmat <- gauss_cor(c(x, xx))
      #kx <- gauss_cor_mat(self$X) + diag(self$nug, self$N)
      kxx <- self$corr_func(XX, theta=self$theta)
      kx.xx <- self$corr_func(self$X, XX, theta=self$theta)

      #mn <- self$pred_mean(XX, kx.xx=kx.xx)
      mn <- pred_meanC(XX, kx.xx, self$mu_hat, self$Kinv, self$Z)
      if (!se.fit & !covmat) {
        return(mn)
      }
      #s2 <- self$pred_var(XX, kxx=kxx, kx.xx=kx.xx)
      #se <- rep(0, length(mn)) # NEG VARS will be 0 for se, NOT SURE I WANT THIS
      #se[s2>=0] <- sqrt(s2[s2>=0])
      if (covmat) {
        #covmatdat <- self$pred_var(XX, kxx=kxx, kx.xx=kx.xx, covmat=T)
        covmatdat <- pred_cov(XX, kxx, kx.xx, self$s2_hat, self$Kinv, self$Z)
        s2 <- diag(covmatdat)
        se <- rep(1e-8, length(mn)) # NEG VARS will be 0 for se, NOT SURE I WANT THIS
        se[s2>=0] <- sqrt(s2[s2>=0])
        return(list(mean=mn, s2=s2, se=se, cov=covmatdat))
      }

      #s2 <- self$pred_var(XX, kxx=kxx, kx.xx=kx.xx, covmat=F)
      s2 <- pred_var(XX, kxx, kx.xx, self$s2_hat, self$Kinv, self$Z)
      se <- rep(0, length(mn)) # NEG VARS will be 0 for se, NOT SURE I WANT THIS
      se[s2>=0] <- sqrt(s2[s2>=0])

      # se.fit but not covmat
      data.frame(mean=mn, s2=s2, se=se)
    },
    pred_mean = function(XX, kx.xx) { # 2-8x faster to use pred_meanC
      c(self$mu_hat + t(kx.xx) %*% self$Kinv %*% (self$Z - self$mu_hat))
    },
    pred_meanC = function(XX, kx.xx) { # Don't use if R uses pass by copy(?)
      pred_meanC(XX, kx.xx, self$mu_hat, self$Kinv, self$Z)
    },
    pred_var = function(XX, kxx, kx.xx, covmat=F) { # 2-4x faster to use C functions pred_var and pred_cov
      self$s2_hat * diag(kxx - t(kx.xx) %*% self$Kinv %*% kx.xx)
    },
    deviance_theta = function (theta) {
      self$deviance(theta=theta, nug=self$nug)
    },
    deviance_theta_log = function (beta) {
      theta <- 10^beta
      self$deviance_theta(theta=theta)
    },
    cool1Dplot = function () {
      if (self$D != 1) stop('Must be 1D')
      minx <- min(self$X)
      maxx <- max(self$X)
      x1 <- minx - .1 * (maxx - minx)
      x2 <- maxx + .1 * (maxx - minx)
      nn <- 201
      x <- seq(x1, x2, length.out = nn)
      px <- self$pred(x, covmat = T)
      n2 <- 20
      Sigma.try <- try(newy <- MASS::mvrnorm(n=n2, mu=px$mean, Sigma=px$cov))
      if (inherits(Sigma.try, "try-error")) {
        message("Adding nugget to cool1Dplot")
        Sigma.try2 <- try(newy <- MASS::mvrnorm(n=n2, mu=px$mean, Sigma=px$cov + diag(self$nug, nrow(px$cov))))
        if (inherits(Sigma.try2, "try-error")) {
          stop("Can't do cool1Dplot")
        }
      }
      plot(x,px$me, type='l', lwd=4, ylim=c(min(newy),max(newy)))
      sapply(1:n2, function(i) points(x, newy[i,], type='l', col='gray'))
      points(self$X, self$Z, pch=19, col=1, cex=2)
    },

    deviance = function (theta=self$theta, nug=self$nug) { #browser()# joint deviance
      K <- self$corr_func(self$X, theta=theta) + diag(nug, self$N)
      Kchol <- try(chol(K))
      if (inherits(Kchol, "try-error")) {return(Inf)}
      Kinv <- chol2inv(Kchol)
      mu_hat <- sum(Kinv %*% self$Z) / sum(Kinv)
      logdetK <- 2 * sum(log(diag(Kchol)))
      logdetK + self$N * log(t(self$Z - mu_hat) %*% (Kinv %*% (self$Z - mu_hat)))
    },
    devianceC2 = function (theta=self$theta, nug=self$nug) { #browser()# joint deviance
      # Not faster than devianceC or even deviance
      K <- self$corr_func(self$X, theta=theta) + diag(nug, self$N)
      Kchol <- try(chol(K))
      if (inherits(Kchol, "try-error")) {return(Inf)}
      Kinv <- chol2inv(Kchol)
      logdetK <- 2 * sum(log(diag(Kchol)))
      logdetK + deviance_part(theta, nug, self$X, self$Z, Kinv)
    },
    devianceC = function (theta=self$theta, nug=self$nug) { #browser()# joint deviance
      # Uses cholC and solveC, not much faster than deviance though, maybe not at all
      # This gives error for n>32, probably from cholC, so DONT USE
      K <- self$corr_func(self$X, theta=theta) + diag(nug, self$N)
      Kchol <- try(cholC(K))
      if (inherits(Kchol, "try-error")) {return(Inf)}
      Kinv <- chol2inv(Kchol)
      mu_hat <- sum(Kinv %*% self$Z) / sum(Kinv)
      logdetK <- 2 * sum(log(diag(Kchol)))
      #logdetK + self$N * log(t(self$Z - mu_hat) %*% solveC(K, self$Z - mu_hat)) # convert 1x1 matrix to num
      logdetK + self$N * log(as.numeric(t(self$Z - mu_hat) %*% (Kinv %*% (self$Z - mu_hat))))
    },
    deviance_log = function (beta=NULL, nug=self$nug, joint=NULL) {#browser()  # joint deviance
      if (!is.null(joint)) {
        beta <- joint[-length(joint)]
        theta <- 10^beta
        nug <- joint[length(joint)]
      } else {
        if (is.null(beta)) theta <- self$theta
        else theta <- 10^beta
      }
      (if (self$useC) self$devianceC2 else self$deviance)(theta=theta, nug=nug)
    },
    deviance_log2 = function (beta=NULL, lognug=NULL, joint=NULL) {#browser()  # joint deviance
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
      (if (self$useC) self$devianceC2 else self$deviance)(theta=theta, nug=nug)
    },
    deviance_grad = function (theta=self$theta, nug=self$nug) { #browser()# joint deviance
      K <- self$corr_func(self$X, theta=theta) + diag(nug, self$N)
      Kchol <- try(chol(K))
      if (inherits(Kchol, "try-error")) {return(Inf)}
      Kinv <- chol2inv(Kchol)
      mu_hat <- sum(Kinv %*% self$Z) / sum(Kinv)
      #logdetK <- 2 * sum(log(diag(Kchol)))
      #logdetK + self$N * log(t(self$Z - mu_hat) %*% (Kinv %*% (self$Z - mu_hat)))

      y <- self$Z - mu_hat
      Kinv.y <- Kinv %*% y
      t2a <- -self$N / (t(y)%*%Kinv.y)
      dD <- rep(NA,self$D)
      for(i in 1:self$D) {
        dK <- K
        for(j in 1:self$N) {
          for(k in 1:self$N) {
            dK[j, k] <- -(self$X[j,i]-self$X[k,i])^2 * dK[j, k]
          }
        }
        #t1 <- sum(diag(Kinv %*% dK))
        t1 <- sum(sapply(1:self$N, function(ii) {sum(Kinv[ii,] * dK[,ii])}))
        t2 <- t2a * t(Kinv.y) %*% dK %*% Kinv.y
        dD[i] <- t1+t2
      }
      dD
    },
    deviance_log_grad = function (beta=NULL, nug=self$nug, joint=NULL) {#browser()  # joint deviance
      if (!is.null(joint)) {
        beta <- joint[-length(joint)]
        theta <- 10^beta
        nug <- joint[length(joint)]
      } else {
        if (is.null(beta)) theta <- self$theta
        else theta <- 10^beta
      }
      self$deviance_grad(theta=theta, nug=nug) * 10^beta * log(10)
    },
    deviance_LLH = function (theta=self$theta, nug=self$nug) { #browser()# joint deviance
      K <- self$corr_func(self$X, theta=theta) + diag(nug, self$N)
      Kchol <- try(chol(K))
      if (inherits(Kchol, "try-error")) {return(Inf)}
      Kinv <- chol2inv(Kchol)
      mu_hat <- sum(Kinv %*% self$Z) / sum(Kinv)
      logdetK <- 2 * sum(log(diag(Kchol)))
      logdetK + (t(self$Z - mu_hat) %*% (Kinv %*% (self$Z - mu_hat)))
    },
    deviance_log_LLH = function (beta=NULL, nug=self$nug, joint=NULL) {#browser()  # joint deviance
      if (!is.null(joint)) {
        beta <- joint[-length(joint)]
        theta <- 10^beta
        nug <- joint[length(joint)]
      } else {
        if (is.null(beta)) theta <- self$theta
        else theta <- 10^beta
      }
      self$deviance_LLH(theta=theta, nug=nug)
    },
    optimOld = function (restarts = 5, theta.update = T, nug.update = self$nug.est) {#browser()
      # Joint MLE search with L-BFGS-B, with restarts
      if (theta.update & nug.update) {
        optim.func <- function(xx) {self$deviance_log(joint=xx)}
      } else if (theta.update & !nug.update) {
        optim.func <- function(xx) {self$deviance_log(beta=xx)}
      } else if (!theta.update & nug.update) {
        optim.func <- function(xx) {self$deviance_log(nug=xx)}
      } else {
        stop("Can't optimize over no variables")
      }
      lower <- c()
      upper <- c()
      start.par <- c()
      start.par0 <- c() # Some default params
      if (theta.update) {
        lower <- c(lower, rep(-5, self$theta_length))
        upper <- c(upper, rep(7, self$theta_length))
        start.par <- c(start.par, log(self$theta, 10))
        start.par0 <- c(start.par0, rep(0, length(self$theta)))
      }
      if (nug.update) {
        lower <- c(lower, self$nug.min)
        upper <- c(upper, Inf)
        start.par <- c(start.par, self$nug)
        start.par0 <- c(start.par0, 1e-6)
      }

      # Find best params with optimization, start with current params in case all give error
      # Current params
      best <- list(par=c(log(self$theta, 10), self$nug), value = self$deviance_log())
      if (self$verbose >= 2) {cat("Optimizing\n");cat("\tInitial values:\n");print(best)}
      details <- data.frame(start=paste(c(self$theta,self$nug),collapse=","),end=NA,value=best$value,func_evals=1,grad_evals=NA,convergence=NA, message=NA, newbest=1)

      # Run optim from current
      current <- try(
        #optim(c(log(self$theta, 10),self$nug), self$deviance_log, method="L-BFGS-B", lower=lower, upper=upper, hessian=F)
        optim(start.par, optim.func, method="L-BFGS-B", lower=lower, upper=upper, hessian=F)
      )
      if (!inherits(current, "try-error")) {#browser()
        #if (self$verbose >= 2) {cat("\tFirst run:\n");print(current)}
        details.new <- data.frame(start=paste(signif(c(self$theta,self$nug),3),collapse=","),end=paste(signif(current$par,3),collapse=","),value=current$value,func_evals=current$counts[1],grad_evals=current$counts[2],convergence=current$convergence, message=current$message, newbest=current$value < best$value, row.names = NULL)
        if (current$value < best$value) {
          best <- current
        }
      } else {#browser() # NEED THESE NEW DETAILS
        if (self$verbose >= 2) {cat("\tFirst run: try-error\n");print(current)}
        details.new <- data.frame(start=paste(signif(c(self$theta,self$nug),3),collapse=","),end="try-error",value=NA,func_evals=NA,grad_evals=NA,convergence=NA, message=current[1], newbest=0)
      }
      details <- rbind(details, details.new)
      if (restarts >= 1) {
        for (i in 1:restarts) {

          if (runif(1) < .33) { # restart near some spot to avoid getting stuck in bad spot
            start.par.i <- start.par0
            #print("start at zero par")
          } else { # jitter from current params
            start.par.i <- start.par
          } #;browser()
          if (theta.update) {start.par.i[1:self$theta_length] <- start.par.i[1:self$theta_length] + rnorm(self$theta_length,0,2)} # jitter betas
          if (nug.update) {start.par.i[length(start.par.i)] <- start.par.i[length(start.par.i)] + rexp(1,1e4)} # jitter nugget
          if (self$verbose >= 2) {cat("\tRestart",i,": starts pars =",start.par.i,"\n")}
          current <- try(
            optim(start.par.i, optim.func, method="L-BFGS-B", lower=lower, upper=upper, hessian=F)
          )
          if (!inherits(current, "try-error")) {
            details.new <- data.frame(start=paste(signif(start.par.i,3),collapse=","),end=paste(signif(current$par,3),collapse=","),value=current$value,func_evals=current$counts[1],grad_evals=current$counts[2],convergence=current$convergence, message=current$message, newbest=current$value < best$value, row.names = NULL)
            if (current$value < best$value) {
              best <- current
            }
          } else{
            details.new <- data.frame(start=paste(signif(c(self$theta,self$nug),3),collapse=","),end="try-error",value=NA,func_evals=NA,grad_evals=NA,convergence=NA, message=current[1], newbest=0)
          }
          details <- rbind(details, details.new)
        }

      }
      if (self$verbose >= 2) {print(details)}
      best
    },
    optimx = function (restarts = 5, theta.update = T, nug.update = self$nug.est) {#browser()
      # Joint MLE search with L-BFGS-B, with restarts
      if (theta.update & nug.update) {
        optim.func <- function(xx) {self$deviance_log(joint=xx)}
      } else if (theta.update & !nug.update) {
        optim.func <- function(xx) {self$deviance_log(beta=xx)}
      } else if (!theta.update & nug.update) {
        optim.func <- function(xx) {self$deviance_log(nug=xx)}
      } else {
        stop("Can't optimize over no variables")
      }
      lower <- c()
      upper <- c()
      start.par <- c()
      start.par0 <- c() # Some default params
      if (theta.update) {
        lower <- c(lower, rep(-5, self$theta_length))
        upper <- c(upper, rep(7, self$theta_length))
        start.par <- c(start.par, log(self$theta, 10))
        start.par0 <- c(start.par0, rep(0, length(self$theta)))
      }
      if (nug.update) {
        lower <- c(lower, self$nug.min)
        upper <- c(upper, Inf)
        start.par <- c(start.par, self$nug)
        start.par0 <- c(start.par0, 1e-6)
      }

      # Find best params with optimization, start with current params in case all give error
      # Current params
      best <- list(par=c(log(self$theta, 10), self$nug), value = self$deviance_log())
      if (self$verbose >= 2) {cat("Optimizing\n");cat("\tInitial values:\n");print(best)}
      details <- data.frame(start=paste(c(self$theta,self$nug),collapse=","),end=NA,value=best$value,func_evals=1,grad_evals=NA,convergence=NA, message=NA, newbest=1)

      # Run optim from current
      current <- try(
        #optim(c(log(self$theta, 10),self$nug), self$deviance_log, method="L-BFGS-B", lower=lower, upper=upper, hessian=F)
        optimx(start.par, optim.func, method="L-BFGS-B", lower=lower, upper=upper, hessian=F)
      )
      if (!inherits(current, "try-error")) {#browser()
        #if (self$verbose >= 2) {cat("\tFirst run:\n");print(current)}
        current.par <- as.numeric((current[paste0("p",1:length(start.par))]))
        details.new <- data.frame(start=paste(signif(c(self$theta,self$nug),3),collapse=","),end=paste(signif(current.par,3),collapse=","),value=current$value,func_evals=current$fevals,grad_evals=current$gevals,convergence=current$convcode, message=current$convcode, newbest=current$value < best$value, row.names = NULL)
        if (current$value < best$value) {
          best <- current
        }
      } else {#browser() # NEED THESE NEW DETAILS
        if (self$verbose >= 2) {cat("\tFirst run: try-error\n");print(current)}
        details.new <- data.frame(start=paste(signif(c(self$theta,self$nug),3),collapse=","),end="try-error",value=NA,func_evals=NA,grad_evals=NA,convergence=NA, message=current[1], newbest=0)
      }
      details <- rbind(details, details.new)
      if (restarts >= 1) {
        for (i in 1:restarts) {

          if (runif(1) < .33) { # restart near some spot to avoid getting stuck in bad spot
            start.par.i <- start.par0
            #print("start at zero par")
          } else { # jitter from current params
            start.par.i <- start.par
          } #;browser()
          if (theta.update) {start.par.i[1:self$theta_length] <- start.par.i[1:self$theta_length] + rnorm(self$theta_length,0,2)} # jitter betas
          if (nug.update) {start.par.i[length(start.par.i)] <- start.par.i[length(start.par.i)] + rexp(1,1e4)} # jitter nugget
          if (self$verbose >= 2) {cat("\tRestart",i,": starts pars =",start.par.i,"\n")}
          current <- try(
            optimx(start.par.i, optim.func, method="L-BFGS-B", lower=lower, upper=upper, hessian=F)
          )
          if (!inherits(current, "try-error")) {#browser()
            current.par <- as.numeric((current[paste0("p",1:length(start.par))]))
            details.new <- data.frame(start=paste(signif(start.par.i,3),collapse=","),end=paste(signif(current.par,3),collapse=","),value=current$value,func_evals=current$fevals,grad_evals=current$gevals,convergence=current$convcode, message=current$convcode, newbest=current$value < best$value, row.names = NULL)
            if (current$value < best$value) {
              best <- current
            }
          } else{
            details.new <- data.frame(start=paste(signif(c(self$theta,self$nug),3),collapse=","),end="try-error",value=NA,func_evals=NA,grad_evals=NA,convergence=NA, message=current[1], newbest=0)
          }
          details <- rbind(details, details.new)
        }

      }
      if (self$verbose >= 2) {print(details)}
      best
    },
    optim = function (restarts = 5, theta.update = T, nug.update = self$nug.est, parallel=self$parallel, parallel.cores=self$parallel.cores) {#browser()
      # Does parallel
      # Joint MLE search with L-BFGS-B, with restarts
      if (theta.update & nug.update) {
        optim.func <- function(xx) {self$deviance_log(joint=xx)}
      } else if (theta.update & !nug.update) {
        optim.func <- function(xx) {self$deviance_log(beta=xx)}
      } else if (!theta.update & nug.update) {
        optim.func <- function(xx) {self$deviance_log(nug=xx)}
      } else {
        stop("Can't optimize over no variables")
      }
      lower <- c()
      upper <- c()
      start.par <- c()
      start.par0 <- c() # Some default params
      if (theta.update) {
        lower <- c(lower, rep(-5, self$theta_length))
        upper <- c(upper, rep(7, self$theta_length))
        start.par <- c(start.par, log(self$theta, 10))
        start.par0 <- c(start.par0, rep(0, length(self$theta)))
      }
      if (nug.update) {
        lower <- c(lower, self$nug.min)
        upper <- c(upper, Inf)
        start.par <- c(start.par, self$nug)
        start.par0 <- c(start.par0, 1e-6)
      }

      # Find best params with optimization, start with current params in case all give error
      # Current params
      best <- list(par=c(log(self$theta, 10), self$nug), value = self$deviance_log())
      if (self$verbose >= 2) {cat("Optimizing\n");cat("\tInitial values:\n");print(best)}
      details <- data.frame(start=paste(c(self$theta,self$nug),collapse=","),end=NA,value=best$value,func_evals=1,grad_evals=NA,convergence=NA, message=NA, stringsAsFactors=F)

      # Run optim from current
      current <- try(
        #optim(c(log(self$theta, 10),self$nug), self$deviance_log, method="L-BFGS-B", lower=lower, upper=upper, hessian=F)
        optim(start.par, optim.func, method="L-BFGS-B", lower=lower, upper=upper, hessian=F)
      )
      if (!inherits(current, "try-error")) {#browser()
        #if (self$verbose >= 2) {cat("\tFirst run:\n");print(current)}
        details.new <- data.frame(start=paste(signif(c(self$theta,self$nug),3),collapse=","),end=paste(signif(current$par,3),collapse=","),value=current$value,func_evals=current$counts[1],grad_evals=current$counts[2],convergence=current$convergence, message=current$message, row.names = NULL, stringsAsFactors=F)
        if (current$value < best$value) {
          best <- current
        }
      } else {#browser() # NEED THESE NEW DETAILS
        if (self$verbose >= 2) {cat("\tFirst run: try-error\n");print(current)}
        details.new <- data.frame(start=paste(signif(c(self$theta,self$nug),3),collapse=","),end="try-error",value=NA,func_evals=NA,grad_evals=NA,convergence=NA, message=current[1], stringsAsFactors=F)

      }
      details <- rbind(details, details.new)
      #if (restarts >= 1) {
        #for (i in 1:restarts) {

        #  details <- rbind(details, details.new)
        #}
        #browser()
      #  restarts.out <- parallel::mclapply(1:restarts, function(i){self$optimRestart(start.par=start.par, start.par0=start.par0, theta.update=theta.update, nug.update=nug.update, optim.func=optim.func, lower=lower, upper=upper)}, mc.cores = self$parallel.cores)
      #  new.details <- t(sapply(restarts.out,function(dd){dd$deta}))
      #  bestparallel <- which.min(sapply(restarts.out,function(i){i$current$val})) #which.min(new.details$value)
      #  if (restarts.out[[bestparallel]]$current$val < best$val) {
      #    best <- restarts.out[[bestparallel]]$current
      #  }
      #  details <- rbind(details, new.details)
      #}
      #browser()
      # runs them in parallel, first starts from current, rest are jittered or random
      if (sys_name == "Windows") {
        # Trying this so it works on Windows
        restarts.out <- lapply( 1:(1+restarts), function(i){self$optimRestart2(start.par=start.par, start.par0=start.par0, theta.update=theta.update, nug.update=nug.update, optim.func=optim.func, lower=lower, upper=upper, jit=(i!=1))})#, mc.cores = parallel.cores)
      } else { # Mac/Unix

        restarts.out <- parallel::mclapply(1:(1+restarts), function(i){self$optimRestart2(start.par=start.par, start.par0=start.par0, theta.update=theta.update, nug.update=nug.update, optim.func=optim.func, lower=lower, upper=upper, jit=(i!=1))}, mc.cores = parallel.cores)
      }
      #restarts.out <- lapply(1:(1+restarts), function(i){self$optimRestart(start.par=start.par, start.par0=start.par0, theta.update=theta.update, nug.update=nug.update, optim.func=optim.func, lower=lower, upper=upper, jit=(i!=1))})
      new.details <- t(sapply(restarts.out,function(dd){dd$deta}))
      bestparallel <- which.min(sapply(restarts.out,function(i){i$current$val})) #which.min(new.details$value)
      if (restarts.out[[bestparallel]]$current$val < best$val) {
        best <- restarts.out[[bestparallel]]$current
      }
      details <- rbind(details, new.details)

      if (self$verbose >= 2) {print(details)}
      best
    },
    optimRestart = function (start.par, start.par0, theta.update, nug.update, optim.func, lower, upper, jit=T) {

      if (runif(1) < .33) { # restart near some spot to avoid getting stuck in bad spot
        start.par.i <- start.par0
        #print("start at zero par")
      } else { # jitter from current params
        start.par.i <- start.par
      } #;browser()
      if (jit) {
        if (theta.update) {start.par.i[1:self$theta_length] <- start.par.i[1:self$theta_length] + rnorm(self$theta_length,0,2)} # jitter betas
        if (nug.update) {start.par.i[length(start.par.i)] <- start.par.i[length(start.par.i)] + rexp(1,1e4)} # jitter nugget
      }
      if (self$verbose >= 2) {cat("\tRestart (parallel): starts pars =",start.par.i,"\n")}
      current <- try(
        optim(start.par.i, optim.func, method="L-BFGS-B", lower=lower, upper=upper, hessian=F)
      )
      if (!inherits(current, "try-error")) {
        details.new <- data.frame(start=paste(signif(start.par.i,3),collapse=","),end=paste(signif(current$par,3),collapse=","),value=current$value,func_evals=current$counts[1],grad_evals=current$counts[2],convergence=current$convergence, message=current$message, row.names = NULL, stringsAsFactors=F)
        #if (current$value < best$value) {
        #  best <- current
        #}
      } else{
        details.new <- data.frame(start=paste(signif(start.par.i,3),collapse=","),end="try-error",value=NA,func_evals=NA,grad_evals=NA,convergence=NA, message=current[1], stringsAsFactors=F)
      }
      list(current=current, details=details.new)
    },
    optim2 = function (restarts = 5, theta.update = T, nug.update = self$nug.est, parallel=self$parallel, parallel.cores=self$parallel.cores) {#browser()
      # Does parallel
      # Joint MLE search with L-BFGS-B, with restarts
      #browser()
      if (theta.update & nug.update) {
        optim.func <- function(xx) {self$deviance_log2(joint=xx)}
      } else if (theta.update & !nug.update) {
        optim.func <- function(xx) {self$deviance_log2(beta=xx)}
      } else if (!theta.update & nug.update) {
        optim.func <- function(xx) {self$deviance_log(lognug=xx)}
      } else {
        stop("Can't optimize over no variables")
      }
      lower <- c()
      upper <- c()
      start.par <- c()
      start.par0 <- c() # Some default params
      if (theta.update) {
        lower <- c(lower, rep(-5, self$theta_length))
        upper <- c(upper, rep(7, self$theta_length))
        start.par <- c(start.par, log(self$theta, 10))
        start.par0 <- c(start.par0, rep(0, length(self$theta)))
      }
      if (nug.update) {
        lower <- c(lower, log(self$nug.min,10))
        upper <- c(upper, Inf)
        start.par <- c(start.par, log(self$nug,10))
        start.par0 <- c(start.par0, -6)
      }

      # Find best params with optimization, start with current params in case all give error
      # Current params
      best <- list(par=c(log(self$theta, 10), log(self$nug,10)), value = self$deviance_log2())
      if (self$verbose >= 2) {cat("Optimizing\n");cat("\tInitial values:\n");print(best)}
      details <- data.frame(start=paste(c(self$theta,self$nug),collapse=","),end=NA,value=best$value,func_evals=1,grad_evals=NA,convergence=NA, message=NA, stringsAsFactors=F)

      # Run optim from current
      current <- try(
        #optim(c(log(self$theta, 10),self$nug), self$deviance_log, method="L-BFGS-B", lower=lower, upper=upper, hessian=F)
        optim(start.par, optim.func, method="L-BFGS-B", lower=lower, upper=upper, hessian=F)
      )
      if (!inherits(current, "try-error")) {#browser()
        #if (self$verbose >= 2) {cat("\tFirst run:\n");print(current)}
        details.new <- data.frame(start=paste(signif(c(self$theta,self$nug),3),collapse=","),end=paste(signif(current$par,3),collapse=","),value=current$value,func_evals=current$counts[1],grad_evals=current$counts[2],convergence=current$convergence, message=current$message, row.names = NULL, stringsAsFactors=F)
        if (current$value < best$value) {
          best <- current
        }
      } else {#browser() # NEED THESE NEW DETAILS
        if (self$verbose >= 2) {cat("\tFirst run: try-error\n");print(current)}
        details.new <- data.frame(start=paste(signif(c(self$theta,self$nug),3),collapse=","),end="try-error",value=NA,func_evals=NA,grad_evals=NA,convergence=NA, message=current[1], stringsAsFactors=F)

      }
      details <- rbind(details, details.new)

      # runs them in parallel, first starts from current, rest are jittered or random
      sys_name <- Sys.info()["sysname"]
      if (sys_name == "Windows") {
        # Trying this so it works on Windows
        restarts.out <- lapply( 1:(1+restarts), function(i){self$optimRestart2(start.par=start.par, start.par0=start.par0, theta.update=theta.update, nug.update=nug.update, optim.func=optim.func, lower=lower, upper=upper, jit=(i!=1))})#, mc.cores = parallel.cores)
      } else { # Mac/Unix
        restarts.out <- parallel::mclapply(1:(1+restarts), function(i){self$optimRestart2(start.par=start.par, start.par0=start.par0, theta.update=theta.update, nug.update=nug.update, optim.func=optim.func, lower=lower, upper=upper, jit=(i!=1))}, mc.cores = parallel.cores)
      }
      new.details <- t(sapply(restarts.out,function(dd){dd$deta}))
      bestparallel <- which.min(sapply(restarts.out,function(i){i$current$val})) #which.min(new.details$value)
      if (restarts.out[[bestparallel]]$current$val < best$val) {
        best <- restarts.out[[bestparallel]]$current
      }
      details <- rbind(details, new.details)

      if (self$verbose >= 2) {print(details)}
      if (nug.update) best$par[length(best$par)] <- 10 ^ (best$par[length(best$par)])
      best
    },
    optimRestart2 = function (start.par, start.par0, theta.update, nug.update, optim.func, lower, upper, jit=T) {
      # FOR lognug RIGHT NOW, seems to be at least as fast, up to 5x on big data, many fewer func_evals
      #    still want to check if it is better or not
      #browser()
      if (runif(1) < .33) { # restart near some spot to avoid getting stuck in bad spot
        start.par.i <- start.par0
        #print("start at zero par")
      } else { # jitter from current params
        start.par.i <- start.par
      } #;browser()
      if (jit) {
        if (theta.update) {start.par.i[1:self$theta_length] <- start.par.i[1:self$theta_length] + rnorm(self$theta_length,0,2)} # jitter betas
        if (nug.update) {start.par.i[length(start.par.i)] <- start.par.i[length(start.par.i)] + min(4, rexp(1,1))} # jitter nugget
      }
      if (self$verbose >= 2) {cat("\tRestart (parallel): starts pars =",start.par.i,"\n")}
      current <- try(
        optim(start.par.i, optim.func, method="L-BFGS-B", lower=lower, upper=upper, hessian=F)
      )
      if (!inherits(current, "try-error")) {
        details.new <- data.frame(start=paste(signif(start.par.i,3),collapse=","),end=paste(signif(current$par,3),collapse=","),value=current$value,func_evals=current$counts[1],grad_evals=current$counts[2],convergence=current$convergence, message=current$message, row.names = NULL, stringsAsFactors=F)
        #if (current$value < best$value) {
        #  best <- current
        #}
      } else{
        details.new <- data.frame(start=paste(signif(start.par.i,3),collapse=","),end="try-error",value=NA,func_evals=NA,grad_evals=NA,convergence=NA, message=current[1], stringsAsFactors=F)
      }
      list(current=current, details=details.new)
    },
    optimBayes = function(theta.update = T, nug.update = F & self$nug.est) {#browser()
      lower <- c()
      upper <- c()
      start.par <- c()
      start.par0 <- c() # Some default params
      if (theta.update) {
        lower <- c(lower, rep(-5, self$theta_length))
        upper <- c(upper, rep(7, self$theta_length))
        start.par <- c(start.par, log(self$theta, 10))
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
      #browser()
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
        if (F) {browser()
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
      #browser()
      bestu <- u[which.min(v),]
      optim(bestu, bgp$pred, lower=-5,upper=7,method="L-BFGS-B")$par
    },
    update = function (Xnew=NULL, Znew=NULL, Xall=NULL, Zall=NULL,
                       restarts = 5, useOptim2=self$useOptim2,
                       theta.update = T, nug.update = self$nug.est) {
      self$update_data(Xnew=Xnew, Znew=Znew, Xall=Xall, Zall=Zall) # Doesn't update Kinv, etc

      pars <- (if(useOptim2) self$optim2 else self$optim)(restarts = restarts, theta.update = theta.update, nug.update = nug.update)$par
      self$nug <- pars[length(pars)]
      self$theta <- 10 ^ pars[-length(pars)]
      self$update_params()

      invisible(self)
    },
    update_data = function(Xnew=NULL, Znew=NULL, Xall=NULL, Zall=NULL) {
      if (!is.null(Xall)) {self$X <- Xall;self$N <- nrow(self$X)} else if (!is.null(Xnew)) {self$X <- rbind(self$X, Xnew);self$N <- nrow(self$X)}
      if (!is.null(Zall)) {self$Z <- Zall} else if (!is.null(Znew)) {self$Z <- rbind(self$Z, Znew)}
      #if (!is.null(Xall) | !is.null(Xnew)) {self$update_params()} # update Kinv, etc, DONT THINK I NEED IT
    },
    update_theta = function (...) {
      self$update(nug.update = F, ...=...)
    },
    update_nugget = function (...) {
      self$update(theta.update = F, ...=...)
    },
    deviance_searchnug = function() {
      optim(self$nug, function(nnug) {self$deviance(theta=self$theta, nug=nnug)}, method="L-BFGS-B", lower=0, upper=Inf, hessian=F)$par
    },
    nugget_update = function () {
      nug <- self$deviance_searchnug()
      #self$theta <- 10 ^ self$deviance_search()$minimum
      self$nug <- nug
      self$update_params()
    },
    grad = function (XX) {#browser()
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
    grad_norm = function (XX) {#browser()
      grad1 <- self$grad(XX)
      if (!is.matrix(grad1)) return(abs(grad1))
      apply(grad1,1, function(xx) {sqrt(sum(xx^2))})
    }#,
    #grad_num = function (XX) {#browser() # NUMERICAL GRAD IS OVER 10 TIMES SLOWER
    #  if (!is.matrix(XX)) {
    #    if (self$D == 1) XX <- matrix(XX, ncol=1)
    #    else if (length(XX) == self$D) XX <- matrix(XX, nrow=1)
    #    else stop('Predict input should be matrix')
    #  } else {
    #    if (ncol(XX) != self$D) {stop("Wrong dimension input")}
    #  }
    #  grad.func <- function(xx) self$pred(xx)$mean
    #  grad.apply.func <- function(xx) numDeriv::grad(grad.func, xx)
    #  grad1 <- apply(XX, 1, grad.apply.func)
    #  if (self$D == 1) return(grad1)
    #  t(grad1)
    #},
    #grad_num_norm = function (XX) {#browser()
    #  grad1 <- self$grad_num(XX)
    #  if (!is.matrix(grad1)) return(abs(grad1))
    #  apply(grad1,1, function(xx) {sqrt(sum(xx^2))})
    #}
    ,
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
