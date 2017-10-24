#' Class providing object with methods for fitting a GP model
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
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
#' @field corr Type of correlation function
#' @field nug.min Minimum value of nugget
#' @field nug Value of the nugget, is estimated unless told otherwise
#' @field theta Length-scale parameters
#' @field separable Are the dimensions separable?
#' @field verbose 0 means nothing printed, 1 prints some, 2 prints most.
#' @field useGrad Should grad be used?
#' @field useC Should C code be used?
#' @field parallel Should the code be run in parallel?
#' @field parallel.cores How many cores are there? It will self detect, do not set yourself.
#' @section Methods:
#' \describe{
#'   \item{\code{new(X, Z, corr="Gauss", verbose=0, separable=T, useC=F,useGrad=T,
#'          parallel=T, useOptim2=T, nug.est=T, ...)}}{This method is used to create object of this class with \code{X} and \code{Z} as the data.}
#'
#'   \item{\code{update(Xnew=NULL, Znew=NULL, Xall=NULL, Zall=NULL,
#' restarts = 5, useOptim2=self$useOptim2,
#' theta.update = T, nug.update = self$nug.est)}}{This method updates the model, adding new data if given, then running optimization again.}
#'   }
GauPro_old <- R6::R6Class(classname = "GauPro_old",
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
    theta = NULL, # theta of length D, will repeat values if nonseparable
    theta_length = NULL, # number of unique theta's
    theta_map = NULL, # Vector of length D, ith component is which theta group ith dimension belongs to, eg c(1,2,1) means first and third dimensions share a theta
    theta_short = NULL, # vector of unique theta's, length D if separable, length 1 if nonsep. Use this for optimization, then transform to full theta using map
    #theta_short_length = NULL,
    separable = NULL,
    mu_hat = NULL,
    s2_hat = NULL,
    K = NULL,
    Kchol = NULL,
    Kinv = NULL,
    verbose = 0,
    useC = TRUE,
    useGrad = FALSE,
    parallel = FALSE,
    parallel.cores = NULL,
    useOptim2 = FALSE,
    deviance_out = NULL, #(theta, nug)
    deviance_grad_out = NULL, #(theta, nug, overwhat)
    deviance_fngr_out = NULL,
    initialize = function(X, Z, corr="Gauss", verbose=0, separable=T, useC=F,useGrad=T,
                          parallel=T, useOptim2=T, nug.est=T,
                          theta_map = NULL,
                          ...) {
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

      if (self$corr == "Gauss") {
        if (self$useC) {
          self$corr_func <- corr_gauss_matrix
          self$deviance_out <- Gaussian_devianceC
          self$deviance_grad_out <- Gaussian_deviance_gradC
          self$deviance_fngr_out <- Gaussian_deviance_fngrC
        } else {
          self$corr_func <- corr_gauss_matrix_noC
          self$deviance_out <- Gaussian_devianceR
          self$deviance_grad_out <- Gaussian_deviance_gradR
          self$deviance_fngr_out <- Gaussian_deviance_fngrR
        }
      } else {
        stop("corr not specified or recognized")
      }
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
    fit = function(X, Z) {
      self$update()
    },
    update_params = function () {
      while(T) {
        self$K <- self$corr_func(self$X, theta=self$theta) + diag(self$nug, self$N)
        try.chol <- try(self$Kchol <- chol(self$K), silent = T)
        if (!inherits(try.chol, "try-error")) {break}
        warning("Can't Cholesky, increasing nugget #7819553")
        oldnug <- self$nug
        self$nug <- max(1e-8, 2 * self$nug)
        print(c(oldnug, self$nug))
      }
      self$Kinv <- chol2inv(self$Kchol)
      self$mu_hat <- sum(self$Kinv %*% self$Z) / sum(self$Kinv)
      self$s2_hat <- c(t(self$Z - self$mu_hat) %*% self$Kinv %*% (self$Z - self$mu_hat) / self$N)
    },
    predict = function(XX, se.fit=F, covmat=F) {
      self$pred(XX=XX, se.fit=se.fit, covmat=covmat)
    },
    pred = function(XX, se.fit=F, covmat=F) {
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
    deviance = function(theta, nug) {
      self$deviance_out(theta, nug, self$X, self$Z)
    },
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
    optim = function (restarts = 5, theta.update = T, nug.update = self$nug.est, parallel=self$parallel, parallel.cores=self$parallel.cores) {
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
        start.par <- c(start.par, log(self$theta_short, 10))
        start.par0 <- c(start.par0, rep(0, self$theta_length))
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

      # runs them in parallel, first starts from current, rest are jittered or random
      sys_name <- Sys.info()["sysname"]
      if (sys_name == "Windows") {
        # Trying this so it works on Windows
        restarts.out <- lapply( 1:(1+restarts), function(i){self$optimRestart(start.par=start.par, start.par0=start.par0, theta.update=theta.update, nug.update=nug.update, optim.func=optim.func, lower=lower, upper=upper, jit=(i!=1))})#, mc.cores = parallel.cores)
      } else { # Mac/Unix

        restarts.out <- parallel::mclapply(1:(1+restarts), function(i){self$optimRestart(start.par=start.par, start.par0=start.par0, theta.update=theta.update, nug.update=nug.update, optim.func=optim.func, lower=lower, upper=upper, jit=(i!=1))}, mc.cores = parallel.cores)
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
    optimRestart = function (start.par, start.par0, theta.update, nug.update, optim.func, lower, upper, jit=T, startAt.par0=F) {

      if (runif(1) < .33 & jit) { # restart near some spot to avoid getting stuck in bad spot
        start.par.i <- start.par0
        #print("start at zero par")
      } else { # jitter from current params
        start.par.i <- start.par
      }
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
    optim2 = function (restarts = 5, theta.update = T, nug.update = self$nug.est, parallel=self$parallel, parallel.cores=self$parallel.cores) {
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
        self$update_params() # This will increase the nugget
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
                       else ii$current$val
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
    optimRestart2 = function (start.par, start.par0, theta.update, nug.update, optim.func, grad.func, optim.fngr, lower, upper, jit=T) {
      # FOR lognug RIGHT NOW, seems to be at least as fast, up to 5x on big data, many fewer func_evals
      #    still want to check if it is better or not
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
      )
      if (self$useGrad) {current$counts <- c(NA,NA);if(is.null(current$message))current$message=NA}
      if (!inherits(current, "try-error")) {
        details.new <- data.frame(start=paste(signif(start.par.i,3),collapse=","),end=paste(signif(current$par,3),collapse=","),value=current$value,func_evals=current$counts[1],grad_evals=current$counts[2],convergence=current$convergence, message=current$message, row.names = NULL, stringsAsFactors=F)
      } else{
        details.new <- data.frame(start=paste(signif(start.par.i,3),collapse=","),end="try-error",value=NA,func_evals=NA,grad_evals=NA,convergence=NA, message=current[1], stringsAsFactors=F)
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
    update = function (Xnew=NULL, Znew=NULL, Xall=NULL, Zall=NULL,
                       restarts = 5, useOptim2=self$useOptim2,
                       theta.update = T, nug.update = self$nug.est) {
      self$update_data(Xnew=Xnew, Znew=Znew, Xall=Xall, Zall=Zall) # Doesn't update Kinv, etc

      pars <- (if(useOptim2) self$optim2 else self$optim)(restarts = restarts, theta.update = theta.update, nug.update = nug.update)$par
      if (nug.update) {self$nug <- pars[length(pars)]}
      if (theta.update) {
        self$theta_short <- 10 ^ pars[1:self$theta_length]
        self$theta <- self$theta_short[self$theta_map]
      }
      self$update_params()

      invisible(self)
    },
    update_data = function(Xnew=NULL, Znew=NULL, Xall=NULL, Zall=NULL) {
      if (!is.null(Xall)) {
        self$X <- if (is.matrix(Xall)) Xall else matrix(Xall,nrow=1)
        self$N <- nrow(self$X)
      } else if (!is.null(Xnew)) {
        self$X <- rbind(self$X, if (is.matrix(Xnew)) Xnew else matrix(Xnew,nrow=1))
        self$N <- nrow(self$X)
      }
      if (!is.null(Zall)) {
        self$Z <- if (is.matrix(Zall))Zall else matrix(Zall,ncol=1)
      } else if (!is.null(Znew)) {
        self$Z <- rbind(self$Z, if (is.matrix(Znew)) Znew else matrix(Znew,ncol=1))
      }
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
    }#,
    #grad_num = function (XX) { # NUMERICAL GRAD IS OVER 10 TIMES SLOWER
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
    #grad_num_norm = function (XX) {
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
