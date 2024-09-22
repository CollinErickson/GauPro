#' Class providing object with methods for fitting a GP model
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @useDynLib GauPro
#' @importFrom Rcpp evalCpp
#' @importFrom stats optim
# @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link[R6]{R6Class}} with methods for fitting GP model.
#' @format \code{\link[R6]{R6Class}} object.
#' @examples
#' #n <- 12
#' #x <- matrix(seq(0,1,length.out = n), ncol=1)
#' #y <- sin(2*pi*x) + rnorm(n,0,1e-1)
#' #gp <- GauPro(X=x, Z=y, parallel=FALSE)
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
#' @field nug.min Minimum value of nugget
#' @field nug Value of the nugget, is estimated unless told otherwise
#' @field verbose 0 means nothing printed, 1 prints some, 2 prints most.
#' @field useGrad Should grad be used?
#' @field useC Should C code be used?
#' @field parallel Should the code be run in parallel?
#' @field parallel_cores How many cores are there? It will self detect,
#' do not set yourself.
#' @field nug.est Should the nugget be estimated?
#' @field param.est Should the parameters be estimated?
#' @field mu_hat Mean estimate
#' @field s2_hat Variance estimate
#' @field K Covariance matrix
#' @field Kchol Cholesky factorization of K
#' @field Kinv Inverse of K
#' @importFrom stats runif
#' @section Methods:
#' \describe{
#'   \item{\code{new(X, Z, corr="Gauss", verbose=0, separable=T, useC=F,useGrad=T,
#'          parallel=T, nug.est=T, ...)}}{This method is used to create object of this class with \code{X} and \code{Z} as the data.}
#'
#'   \item{\code{update(Xnew=NULL, Znew=NULL, Xall=NULL, Zall=NULL,
#' restarts = 5,
#' param_update = T, nug.update = self$nug.est)}}{This method updates the model, adding new data if given, then running optimization again.}
#'   }
GauPro_base <- R6::R6Class(
  classname = "GauPro",
  public = list(
    X = NULL,
    Z = NULL,
    N = NULL,
    D = NULL,
    nug = NULL,
    nug.min = NULL,
    nug.est = NULL,
    param.est = NULL, # Whether parameters besides nugget (theta) should be updated
    mu_hat = NULL,
    s2_hat = NULL,
    #' @description Correlation function
    #' @param ... Does nothing
    corr_func = function(...){}, # When this was NULL the child didn't overwrite with own method, it stayed as NULL
    K = NULL,
    Kchol = NULL,
    Kinv = NULL,
    verbose = 0,
    useC = TRUE,
    useGrad = FALSE,
    parallel = NULL,
    parallel_cores = NULL,
    #deviance_out = NULL, #(theta, nug)
    #deviance_grad_out = NULL, #(theta, nug, overwhat)
    #deviance_fngr_out = NULL,
    #' @description Create GauPro object
    #' @param X Matrix whose rows are the input points
    #' @param Z Output points corresponding to X
    #' @param verbose Amount of stuff to print. 0 is little, 2 is a lot.
    #' @param useC Should C code be used when possible? Should be faster.
    #' @param useGrad Should the gradient be used?
    #' @param parallel Should code be run in parallel? Make optimization
    #' faster but uses more computer resources.
    #' @param nug Value for the nugget. The starting value if estimating it.
    #' @param nug.min Minimum allowable value for the nugget.
    #' @param nug.est Should the nugget be estimated?
    #' @param param.est Should the kernel parameters be estimated?
    #' @param ... Not used
    initialize = function(X, Z, verbose=0, useC=F,useGrad=T,
                          parallel=FALSE,
                          nug=1e-6, nug.min=1e-8, nug.est=T,
                          param.est = TRUE,
                          ...) {
      #self$initialize_GauPr(X=X,Z=Z,verbose=verbose,useC=useC,useGrad=useGrad,
      #                      parallel=parallel, nug.est=nug.est)
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

      self$nug <- nug
      self$nug.min <- nug.min
      self$nug.est <- nug.est
      self$param.est <- param.est
      self$useC <- useC
      self$useGrad <- useGrad
      self$parallel <- parallel
      if (self$parallel) {self$parallel_cores <- parallel::detectCores()}
      else {self$parallel_cores <- 1}


      invisible(self)
    },
    #' @description Not used
    initialize_GauPr = function() {
    },
    #' @description Fit the model, never use this function
    #' @param X Not used
    #' @param Z Not used
    fit = function(X, Z) {
      self$update()
    },
    #' @description Update Covariance matrix and estimated parameters
    update_K_and_estimates = function () {
      # Update K, Kinv, mu_hat, and s2_hat, maybe nugget too
      while(T) {
        self$K <- self$corr_func(self$X) + diag(self$nug, self$N)
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
    #' @description Predict mean and se for given matrix
    #' @param XX Points to predict at
    #' @param se.fit Should the se be returned?
    #' @param covmat Should the covariance matrix be returned?
    #' @param split_speed Should the predictions be split up for speed
    predict = function(XX, se.fit=F, covmat=F, split_speed=T) {
      self$pred(XX=XX, se.fit=se.fit, covmat=covmat, split_speed=split_speed)
    },
    #' @description Predict mean and se for given matrix
    #' @param XX Points to predict at
    #' @param se.fit Should the se be returned?
    #' @param covmat Should the covariance matrix be returned?
    #' @param split_speed Should the predictions be split up for speed
    pred = function(XX, se.fit=F, covmat=F, split_speed=T) {
      if (!is.matrix(XX)) {
        if (self$D == 1) XX <- matrix(XX, ncol=1)
        else if (length(XX) == self$D) XX <- matrix(XX, nrow=1)
        else stop('Predict input should be matrix')
      }

      N <- nrow(XX)
      # Split speed makes predictions for groups of rows separately.
      # Fastest is for about 40.
      if (split_speed & N >= 200 & !covmat) {#print('In split speed')
        mn <- numeric(N)
        if (se.fit) {
          s2 <- numeric(N)
          se <- numeric(N)
          #se <- rep(0, length(mn)) # NEG VARS will be 0 for se, NOT SURE I WANT THIS
        }

        ni <- 40 # batch size
        Nni <- ceiling(N/ni)-1
        for (j in 0:Nni) {
          XXj <- XX[(j*ni+1):(min((j+1)*ni,N)), , drop=FALSE]
          # kxxj <- self$corr_func(XXj)
          # kx.xxj <- self$corr_func(self$X, XXj)
          predj <- self$pred_one_matrix(XX=XXj, se.fit=se.fit, covmat=covmat)
          #mn[(j*ni+1):(min((j+1)*ni,N))] <- pred_meanC(XXj, kx.xxj, self$mu_hat, self$Kinv, self$Z)
          if (!se.fit) { # if no se.fit, just set vector
            mn[(j*ni+1):(min((j+1)*ni,N))] <- predj
          } else { # otherwise set all three from data.frame
            mn[(j*ni+1):(min((j+1)*ni,N))] <- predj$mean
            #s2j <- pred_var(XXj, kxxj, kx.xxj, self$s2_hat, self$Kinv, self$Z)
            #s2[(j*ni+1):(min((j+1)*ni,N))] <- s2j
            s2[(j*ni+1):(min((j+1)*ni,N))] <- predj$s2
            se[(j*ni+1):(min((j+1)*ni,N))] <- predj$se

          }
        }
        #se[s2>=0] <- sqrt(s2[s2>=0])
        if (!se.fit) {# covmat is always FALSE for split_speed } & !covmat) {
          return(mn)
        } else {
          return(data.frame(mean=mn, s2=s2, se=se))
        }

      } else {
        return(self$pred_one_matrix(XX=XX, se.fit=se.fit, covmat=covmat))
      }
    },
    #' @description Predict mean and se for given matrix
    #' @param XX Points to predict at
    #' @param se.fit Should the se be returned?
    #' @param covmat Should the covariance matrix be returned?
    pred_one_matrix = function(XX, se.fit=F, covmat=F) {
      # input should already be check for matrix
      kxx <- self$corr_func(XX) + self$nug
      kx.xx <- self$corr_func(self$X, XX)
      mn <- pred_meanC(XX, kx.xx, self$mu_hat, self$Kinv, self$Z)

      if (!se.fit & !covmat) {
        return(mn)
      }
      if (covmat) {
        #covmatdat <- self$pred_var(XX, kxx=kxx, kx.xx=kx.xx, covmat=T)
        covmatdat <- pred_cov(XX, kxx, kx.xx, self$s2_hat, self$Kinv, self$Z)
        s2 <- diag(covmatdat)
        se <- rep(1e-8, length(mn)) # NEG VARS will be 0 for se, NOT SURE I WANT THIS
        se[s2>=0] <- sqrt(s2[s2>=0])
        return(list(mean=mn, s2=s2, se=se, cov=covmatdat))
      }

      s2 <- pred_var(XX, kxx, kx.xx, self$s2_hat, self$Kinv, self$Z)
      se <- rep(0, length(mn)) # NEG VARS will be 0 for se, NOT SURE I WANT THIS
      se[s2>=0] <- sqrt(s2[s2>=0])

      # se.fit but not covmat
      data.frame(mean=mn, s2=s2, se=se)
    },
    #' @description Predict mean
    #' @param XX Points to predict at
    #' @param kx.xx Covariance matrix between X and XX
    pred_mean = function(XX, kx.xx) { # 2-8x faster to use pred_meanC
      c(self$mu_hat + t(kx.xx) %*% self$Kinv %*% (self$Z - self$mu_hat))
    },
    #' @description Predict mean using C code
    #' @param XX Points to predict at
    #' @param kx.xx Covariance matrix between X and XX
    pred_meanC = function(XX, kx.xx) { # Don't use if R uses pass by copy(?)
      pred_meanC(XX, kx.xx, self$mu_hat, self$Kinv, self$Z)
    },
    #' @description Predict variance
    #' @param XX Points to predict at
    #' @param kxx Covariance matrix of XX with itself
    #' @param kx.xx Covariance matrix between X and XX
    #' @param covmat Not used
    pred_var = function(XX, kxx, kx.xx, covmat=F) { # 2-4x faster to use C functions pred_var and pred_cov
      self$s2_hat * diag(kxx - t(kx.xx) %*% self$Kinv %*% kx.xx)
    },
    #' @description Predict at X using leave-one-out. Can use for diagnostics.
    #' @param se.fit Should the standard error and t values be returned?
    pred_LOO = function(se.fit=FALSE) {
      # Predict LOO (leave-one-out) on data used to fit model
      # See vignette for explanation of equations
      # If se.fit==T, then calculate the LOO se and the corresponding t score
      Z_LOO <- numeric(self$N)
      if (se.fit) {Z_LOO_se <- numeric(self$N)}
      Z_trend <- self$mu_hat #self$trend$Z(self$X)
      for (i in 1:self$N) {
        E <- self$Kinv[-i, -i] # Kinv without i
        b <- self$K[    i, -i] # K    between i and rest
        g <- self$Kinv[ i, -i] # Kinv between i and rest
        Ainv <- E + E %*% b %*% g / (1-sum(g*b)) # Kinv for K if i wasn't in K
        Zi_LOO <- Z_trend + c(b %*% Ainv %*% (self$Z[-i] - Z_trend))
        Z_LOO[i] <- Zi_LOO
        if (se.fit) { # Need to use s2_hat
          Zi_LOO_se <- sqrt(self$K[i,i] - c(b %*% Ainv %*% b))
          Z_LOO_se[i] <- Zi_LOO_se
        }
      }
      if (se.fit) { # Return df with se and t if se.fit
        Z_LOO_se <- Z_LOO_se * sqrt(self$s2_hat)
        t_LOO <- (self$Z - Z_LOO) / Z_LOO_se
        data.frame(fit=Z_LOO, se.fit=Z_LOO_se, t=t_LOO)
      } else { # Else just mean LOO
        Z_LOO
      }
    },
    #' @description Plot the object
    #' @param ... Parameters passed to cool1Dplot(), plot2D(), or plotmarginal()
    plot = function(...) {
      if (self$D == 1) {
        self$cool1Dplot(...)
      } else if (x$D == 2) {
        self$plot2D(...)
      } else {
        stop("No plot method for higher than 2 dimension")
        # self$plotmarginal(...)
      }
    },
    #' @description Make cool 1D plot
    #' @param n2 Number of things to plot
    #' @param nn Number of things to plot
    #' @param col2 color
    #' @param ylab y label
    #' @param xlab x label
    #' @param xmin xmin
    #' @param xmax xmax
    #' @param ymax ymax
    #' @param ymin ymin
    cool1Dplot = function (n2=20, nn=201, col2="gray",
                           xlab='x', ylab='y',
                           xmin=NULL, xmax=NULL,
                           ymin=NULL, ymax=NULL
    ) {
      if (self$D != 1) stop('Must be 1D')
      # Letting user pass in minx and maxx
      if (is.null(xmin)) {
        minx <- min(self$X)
      } else {
        minx <- xmin
      }
      if (is.null(xmax)) {
        maxx <- max(self$X)
      } else {
        maxx <- xmax
      }
      # minx <- min(self$X)
      # maxx <- max(self$X)
      x1 <- minx - .1 * (maxx - minx)
      x2 <- maxx + .1 * (maxx - minx)
      # nn <- 201
      x <- seq(x1, x2, length.out = nn)
      px <- self$pred(x, covmat = T)
      # n2 <- 20
      Sigma.try <- try(newy <- MASS::mvrnorm(n=n2, mu=px$mean, Sigma=px$cov))
      if (inherits(Sigma.try, "try-error")) {
        message("Adding nugget to cool1Dplot")
        Sigma.try2 <- try(newy <- MASS::mvrnorm(n=n2, mu=px$mean, Sigma=px$cov + diag(self$nug, nrow(px$cov))))
        if (inherits(Sigma.try2, "try-error")) {
          stop("Can't do cool1Dplot")
        }
      }
      # plot(x,px$me, type='l', lwd=4, ylim=c(min(newy),max(newy)),
      #      xlab=xlab, ylab=ylab)
      # sapply(1:n2, function(i) points(x, newy[i,], type='l', col=col2))
      # points(self$X, self$Z, pch=19, col=1, cex=2)

      # Setting ylim, giving user option
      if (is.null(ymin)) {
        miny <- min(newy)
      } else {
        miny <- ymin
      }
      if (is.null(ymax)) {
        maxy <- max(newy)
      } else {
        maxy <- ymax
      }

      # Redo to put gray lines on bottom
      for (i in 1:n2) {
        if (i == 1) {
          plot(x, newy[i,], type='l', col=col2,
               # ylim=c(min(newy),max(newy)),
               ylim=c(miny,maxy),
               xlab=xlab, ylab=ylab)
        } else {
          points(x, newy[i,], type='l', col=col2)
        }
      }
      points(x,px$me, type='l', lwd=4)
      points(self$X, self$Z, pch=19, col=1, cex=2)
    },
    #' @description Make 1D plot
    #' @param n2 Number of things to plot
    #' @param nn Number of things to plot
    #' @param col2 Color of the prediction interval
    #' @param ylab y label
    #' @param xlab x label
    #' @param xmin xmin
    #' @param xmax xmax
    #' @param ymax ymax
    #' @param ymin ymin
    plot1D = function(n2=20, nn=201, col2=2, #"gray",
                      xlab='x', ylab='y',
                      xmin=NULL, xmax=NULL,
                      ymin=NULL, ymax=NULL) {
      if (self$D != 1) stop('Must be 1D')
      # Letting user pass in minx and maxx
      if (is.null(xmin)) {
        minx <- min(self$X)
      } else {
        minx <- xmin
      }
      if (is.null(xmax)) {
        maxx <- max(self$X)
      } else {
        maxx <- xmax
      }
      # minx <- min(self$X)
      # maxx <- max(self$X)
      x1 <- minx - .1 * (maxx - minx)
      x2 <- maxx + .1 * (maxx - minx)
      # nn <- 201
      x <- seq(x1, x2, length.out = nn)
      px <- self$pred(x, se=T)
      # n2 <- 20

      # Setting ylim, giving user option
      if (is.null(ymin)) {
        miny <- min(px$mean - 2*px$se)
      } else {
        miny <- ymin
      }
      if (is.null(ymax)) {
        maxy <- max(px$mean + 2*px$se)
      } else {
        maxy <- ymax
      }

      plot(x, px$mean+2*px$se, type='l', col=col2, lwd=2,
           # ylim=c(min(newy),max(newy)),
           ylim=c(miny,maxy),
           xlab=xlab, ylab=ylab)
      points(x, px$mean-2*px$se, type='l', col=col2, lwd=2)
      points(x,px$me, type='l', lwd=4)
      points(self$X,
             # if (self$normalize) {self$Z * self$normalize_sd + self$normalize_mean}
             # else {self$Z},
             self$Z,
             pch=19, col=1, cex=2)
    },
    #' @description Make 2D plot
    plot2D = function() {
      if (self$D != 2) {stop("plot2D only works in 2D")}
      mins <- apply(self$X, 2, min)
      maxs <- apply(self$X, 2, max)
      xmin <- mins[1] - .03 * (maxs[1] - mins[1])
      xmax <- maxs[1] + .03 * (maxs[1] - mins[1])
      ymin <- mins[2] - .03 * (maxs[2] - mins[2])
      ymax <- maxs[2] + .03 * (maxs[2] - mins[2])
      ContourFunctions_cf_func(self$predict, batchmax=Inf,
                                xlim=c(xmin, xmax),
                                ylim=c(ymin, ymax),
                                pts=self$X)
    },
    #' @description Calculate the log likelihood, don't use this
    #' @param mu Mean vector
    #' @param s2 s2 param
    loglikelihood = function(mu=self$mu_hat, s2=self$s2_hat) {
      -.5 * (self$N*log(s2) + log(det(self$K)) +
               t(self$Z - mu)%*%self$Kinv%*%(self$Z - mu)/s2)
    },
    #' @description Optimize parameters
    #' @param restarts Number of restarts to do
    #' @param param_update Should parameters be updated?
    #' @param nug.update Should nugget be updated?
    #' @param parallel Should restarts be done in parallel?
    #' @param parallel_cores If running parallel, how many cores should be used?
    optim = function (restarts = 5, param_update = T, nug.update = self$nug.est,
                      parallel=self$parallel,
                      parallel_cores=self$parallel_cores) {
      # Does parallel
      # Joint MLE search with L-BFGS-B, with restarts
      #if (param_update & nug.update) {
      #  optim.func <- function(xx) {self$deviance_log2(joint=xx)}
      #  grad.func <- function(xx) {self$deviance_log2_grad(joint=xx)}
      #  optim.fngr <- function(xx) {self$deviance_log2_fngr(joint=xx)}
      #} else if (param_update & !nug.update) {
      #  optim.func <- function(xx) {self$deviance_log2(beta=xx)}
      #  grad.func <- function(xx) {self$deviance_log2_grad(beta=xx)}
      #  optim.fngr <- function(xx) {self$deviance_log2_fngr(beta=xx)}
      #} else if (!param_update & nug.update) {
      #  optim.func <- function(xx) {self$deviance_log2(lognug=xx)}
      #  grad.func <- function(xx) {self$deviance_log2_grad(lognug=xx)}
      #  optim.fngr <- function(xx) {self$deviance_log2_fngr(lognug=xx)}
      #} else {
      #  stop("Can't optimize over no variables")
      #}
      optim_functions <- self$get_optim_functions(param_update=param_update, nug.update=nug.update)
      #optim.func <- self$get_optim_func(param_update=param_update, nug.update=nug.update)
      #optim.grad <- self$get_optim_grad(param_update=param_update, nug.update=nug.update)
      #optim.fngr <- self$get_optim_fngr(param_update=param_update, nug.update=nug.update)
      optim.func <- optim_functions[[1]]
      optim.grad <- optim_functions[[2]]
      optim.fngr <- optim_functions[[3]]


      # Set starting parameters and bounds
      lower <- c()
      upper <- c()
      start.par <- c()
      start.par0 <- c() # Some default params
      if (param_update) {
        lower <- c(lower, self$param_optim_lower())#rep(-5, self$theta_length))
        upper <- c(upper, self$param_optim_upper())#rep(7, self$theta_length))
        start.par <- c(start.par, self$param_optim_start())#log(self$theta_short, 10))
        start.par0 <- c(start.par0, self$param_optim_start0())#rep(0, self$theta_length))
      }
      if (nug.update) {
        lower <- c(lower, log(self$nug.min,10))
        upper <- c(upper, Inf)
        start.par <- c(start.par, log(self$nug,10))
        start.par0 <- c(start.par0, -6)
      }


      # This will make sure it at least can start
      # Run before it sets initial parameters
      try.devlog <- try(devlog <- optim.func(start.par), silent = T)
      if (inherits(try.devlog, "try-error")) {
        warning("Current nugget doesn't work, increasing it #31973")
        self$update_K_and_estimates() # This will increase the nugget until cholesky works
        devlog <- optim.func(start.par)
      }

      # Find best params with optimization, start with current params in case all give error
      # Current params
      #best <- list(par=c(log(self$theta_short, 10), log(self$nug,10)), value = devlog)
      best <- list(par=start.par, value = devlog)
      if (self$verbose >= 2) {
        cat("Optimizing\n");cat("\tInitial values:\n");print(best)
      }
      details <- data.frame(start=paste(start.par,collapse=","),end=NA,
                            value=best$value,func_evals=1,grad_evals=NA,
                            convergence=NA, message=NA, stringsAsFactors=F)


      # runs them in parallel, first starts from current, rest are jittered or random
      sys_name <- Sys.info()["sysname"]
      if (sys_name == "Windows" | !self$parallel) {
        # Trying this so it works on Windows
        restarts.out <- lapply(
          1:(1+restarts),
          function(i){
            self$optimRestart(
              start.par=start.par, start.par0=start.par0,
              param_update=param_update, nug.update=nug.update,
              optim.func=optim.func, optim.grad=optim.grad,
              optim.fngr=optim.fngr, lower=lower, upper=upper,
              jit=(i!=1))})#, mc.cores = parallel_cores)
      } else { # Mac/Unix
        restarts.out <- parallel::mclapply(
          1:(1+restarts),
          function(i){
            self$optimRestart(
              start.par=start.par,
              start.par0=start.par0,
              param_update=param_update,
              nug.update=nug.update,
              optim.func=optim.func,
              optim.grad=optim.grad,
              optim.fngr=optim.fngr,
              lower=lower, upper=upper,
              jit=(i!=1))}, mc.cores = parallel_cores)
      }
      new.details <- t(sapply(restarts.out,
                              function(dd){dd$deta}))
      vals <- sapply(restarts.out,
                     function(ii){
                       if (inherits(ii$current,"try-error")){Inf}
                       else ii$current$val
                     }
      )
      bestparallel <- which.min(vals) #which.min(new.details$value)
      if(inherits(try(restarts.out[[bestparallel]]$current$val, silent = T),
                  "try-error")) {
        # need this in case all are restart vals are Inf
        print("All restarts had error, keeping initial")
      } else if (restarts.out[[bestparallel]]$current$val < best$val) {
        best <- restarts.out[[bestparallel]]$current
      }
      details <- rbind(details, new.details)

      if (self$verbose >= 2) {print(details)}

      # If new nug is below nug.min, optimize again with fixed nug
      # Moved into update_params, since I don't want to set nugget here

      if (nug.update) best$par[length(best$par)] <- 10 ^ (best$par[length(best$par)])
      best
    },
    #' @description Run a single optimization restart.
    #' @param start.par Starting parameters
    #' @param start.par0 Starting parameters
    #' @param param_update Should parameters be updated?
    #' @param nug.update Should nugget be updated?
    #' @param optim.func Function to optimize.
    #' @param optim.grad Gradient of function to optimize.
    #' @param optim.fngr Function that returns the function value
    #' and its gradient.
    #' @param lower Lower bounds for optimization
    #' @param upper Upper bounds for optimization
    #' @param jit Is jitter being used?
    optimRestart = function (start.par, start.par0, param_update, nug.update,
                             optim.func, optim.grad, optim.fngr, lower, upper,
                             jit=T) {
      # FOR lognug RIGHT NOW, seems to be at least as fast, up to 5x on big data, many fewer func_evals
      #    still want to check if it is better or not

      if (runif(1) < .33 & jit) { # restart near some spot to avoid getting stuck in bad spot
        start.par.i <- start.par0
        #print("start at zero par")
      } else { # jitter from current params
        start.par.i <- start.par
      }
      if (jit) {
        #if (param_update) {start.par.i[1:self$theta_length] <- start.par.i[1:self$theta_length] +
        # rnorm(self$theta_length,0,2)} # jitter betas
        theta_indices <- 1:length(self$param_optim_start()) #if () -length(start.par.i)
        if (param_update) {
          start.par.i[theta_indices] <- start.par.i[theta_indices] +
            self$param_optim_jitter(start.par.i[theta_indices])
        } # jitter betas
        if (nug.update) {
          start.par.i[length(start.par.i)] <- start.par.i[length(start.par.i)] +
            min(4, rexp(1,1))} # jitter nugget
      }
      if (self$verbose >= 2) {cat("\tRestart (parallel): starts pars =",start.par.i,"\n")}
      current <- try(
        if (self$useGrad) {
          if (is.null(optim.fngr)) {
            lbfgs::lbfgs(optim.func, optim.grad, start.par.i, invisible=1)
          } else {
            lbfgs_share(optim.fngr, start.par.i, invisible=1) # 1.7x speedup uses grad_share
          }
        } else {
          optim(start.par.i, optim.func, method="L-BFGS-B",
                lower=lower, upper=upper, hessian=F)
        }
      )
      if (!inherits(current, "try-error")) {
        if (self$useGrad) {
          current$counts <- c(NA,NA);
          if(is.null(current$message))current$message=NA
        }
        details.new <- data.frame(
          start=paste(signif(start.par.i,3),collapse=","),
          end=paste(signif(current$par,3),collapse=","),
          value=current$value,func_evals=current$counts[1],
          grad_evals=current$counts[2],convergence=current$convergence,
          message=current$message, row.names = NULL, stringsAsFactors=F)
      } else{
        details.new <- data.frame(
          start=paste(signif(start.par.i,3),collapse=","),
          end="try-error",value=NA,func_evals=NA,grad_evals=NA,
          convergence=NA, message=current[1], stringsAsFactors=F)
      }
      list(current=current, details=details.new)
    },
    #' @description Update the model, can be data and parameters
    #' @param Xnew New X matrix
    #' @param Znew New Z values
    #' @param Xall Matrix with all X values
    #' @param Zall All Z values
    #' @param restarts Number of optimization restarts
    #' @param param_update Should the parameters be updated?
    #' @param nug.update Should the nugget be updated?
    #' @param no_update Should none of the parameters/nugget be updated?
    update = function(Xnew=NULL, Znew=NULL, Xall=NULL, Zall=NULL,
                      restarts = 5,
                      param_update = self$param.est, nug.update = self$nug.est,
                      no_update=FALSE) {
      self$update_data(Xnew=Xnew, Znew=Znew, Xall=Xall, Zall=Zall)
      # Doesn't update Kinv, etc

      if (!no_update || (!param_update && !nug.update)) {
        # This option lets it skip parameter optimization entirely
        self$update_params(restarts=restarts, param_update=param_update,
                           nug.update=nug.update)
      }

      self$update_K_and_estimates()

      invisible(self)
    },
    #' @description Update the data
    #' @param Xnew New X matrix
    #' @param Znew New Z values
    #' @param Xall Matrix with all X values
    #' @param Zall All Z values
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
      #if (!is.null(Xall) | !is.null(Xnew)) {self$update_K_and_estimates()} # update Kinv, etc, DONT THINK I NEED IT
    },
    #' @description Update the correlation parameters
    #' @param ... Args passed to update
    update_corrparams = function (...) {
      self$update(nug.update = F, ...=...)
    },
    #' @description Update the nugget
    #' @param ... Args passed to update
    update_nugget = function (...) {
      self$update(param_update = F, ...=...)
    },
    #' @description Optimize deviance for nugget
    deviance_searchnug = function() {
      optim(self$nug, function(nnug) {self$deviance(nug=nnug)},
            method="L-BFGS-B", lower=0, upper=Inf, hessian=F)$par
    },
    #' @description Update the nugget
    nugget_update = function () {
      nug <- self$deviance_searchnug()
      self$nug <- nug
      self$update_K_and_estimates()
    },
    #' @description Calculate the norm of the gradient at XX
    #' @param XX Points to calculate at
    grad_norm = function (XX) {
      grad1 <- self$grad(XX)
      if (!is.matrix(grad1)) return(abs(grad1))
      apply(grad1,1, function(xx) {sqrt(sum(xx^2))})
    },
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
    #},
    #' @description Sample at XX
    #' @param XX Input points to sample at
    #' @param n Number of samples
    sample = function(XX, n=1) {
      # Generates n samples at rows of XX
      px <- self$pred(XX, covmat = T)
      Sigma.try <- try(newy <- MASS::mvrnorm(n=n, mu=px$mean, Sigma=px$cov))
      if (inherits(Sigma.try, "try-error")) {
        message("Adding nugget to get sample")
        Sigma.try2 <- try(newy <- MASS::mvrnorm(n=n, mu=px$mean, Sigma=px$cov + diag(self$nug, nrow(px$cov))))
        if (inherits(Sigma.try2, "try-error")) {
          stop("Can't do sample, can't factor Sigma")
        }
      }
      newy # Not transposing matrix since it gives var a problem
    },
    #' @description Print object
    print = function() {
      cat("GauPro object\n")
      cat(paste0("\tD = ", self$D, ", N = ", self$N,"\n"))
      cat(paste0("\tNugget = ", signif(self$nug, 3), "\n"))
      cat("\tRun update to add data and/or optimize again\n")
      cat("\tUse pred to get predictions at new points\n")
      invisible(self)
    }
  ),
  private = list(

  )
)
