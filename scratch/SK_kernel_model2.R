#' GauPro model that uses kernels
#'
#' Class providing object with methods for fitting a GP model.
#' Allows for different kernel and trend functions to be used.
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
#' gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian$new(1),
#'                               parallel=FALSE)
#' gp$predict(.454)
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
#' @field corr Type of correlation function
#' @field nug.min Minimum value of nugget
#' @field nug Value of the nugget, is estimated unless told otherwise
#' @field separable Are the dimensions separable?
#' @field verbose 0 means nothing printed, 1 prints some, 2 prints most.
#' @field useGrad Should grad be used?
#' @field useC Should C code be used?
#' @field parallel Should the code be run in parallel?
#' @field parallel_cores How many cores are there? By default it detects.
#' @section Methods:
#' \describe{
#'   \item{\code{new(X, Z, corr="Gauss", verbose=0, separable=T, useC=F,
#'                   useGrad=T,
#'          parallel=T, nug.est=T, ...)}}{
#'          This method is used to create object of this
#'          class with \code{X} and \code{Z} as the data.}
#'
#'   \item{\code{update(Xnew=NULL, Znew=NULL, Xall=NULL, Zall=NULL,
#' restarts = 5,
#' param_update = T, nug.update = self$nug.est)}}{This method updates the
#' model, adding new data if given, then running optimization again.}
#'   }
GauProSK_kernel_model <- R6::R6Class(
  classname = "GauProSK",
  public = list(
    X = NULL,
    Z = NULL,
    N = NULL,
    D = NULL,
    kernel = NULL,
    trend = NULL,
    nug = NULL,
    nug.min = NULL,
    nug.max = NULL,
    nug.est = NULL,
    param.est = NULL,
    # Whether parameters besides nugget (theta) should be updated
    # mu_hat = NULL,
    mu_hatX = NULL,
    s2_hat = NULL,
    K = NULL,
    Kchol = NULL,
    Kinv = NULL,
    Kinv_Z_minus_mu_hatX = NULL,
    verbose = 0,
    useC = TRUE,
    useGrad = FALSE,
    parallel = NULL,
    parallel_cores = NULL,
    restarts = NULL,
    normalize = NULL,
    # Should the Z values be normalized for internal computations?
    normalize_mean = NULL,
    normalize_sd = NULL,
    optimizer = NULL, # L-BFGS-B, BFGS
    #deviance_out = NULL, #(theta, nug)
    #deviance_grad_out = NULL, #(theta, nug, overwhat)
    #deviance_fngr_out = NULL,
    initialize = function(X, Z, A,
                          kernel, trend,
                          verbose=0, useC=F,useGrad=T,
                          parallel=FALSE, parallel_cores="detect",
                          nug=1e-6, nug.min=1e-8, nug.max=Inf, nug.est=TRUE,
                          param.est = TRUE, restarts = 5,
                          normalize = FALSE, optimizer="L-BFGS-B",
                          ...) {
      #self$initialize_GauPr(X=X,Z=Z,verbose=verbose,useC=useC,
      #                      useGrad=useGrad,
      #                      parallel=parallel, nug.est=nug.est)
      self$X <- X
      self$Z <- matrix(Z, ncol=1)
      self$A <- A
      self$normalize <- normalize
      if (self$normalize) {
        self$normalize_mean <- mean(self$Z)
        self$normalize_sd <- sd(self$Z)
        self$Z <- (self$Z - self$normalize_mean) / self$normalize_sd
      }
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

      # Set kernel
      if ("R6ClassGenerator" %in% class(kernel)) {
        # Let generator be given so D can be set auto
        self$kernel <- kernel$new(D=self$D)
      } else if ("GauPro_kernel" %in% class(kernel)) {
        # Otherwise it should already be a kernel
        self$kernel <- kernel
      } else {
        stop("Error: bad kernel #68347")
      }

      # Set trend
      if (missing(trend)) {
        self$trend <- trend_c$new()
      } else if ("GauPro_trend" %in% class(trend)) {
        self$trend <- trend
      } else if ("R6ClassGenerator" %in% class(trend)) {
        self$trend <- trend$new(D=self$D)
      }

      self$nug <- min(max(nug, nug.min), nug.max)
      self$nug.min <- nug.min
      self$nug.max <- nug.max
      self$nug.est <- nug.est
      # if (nug.est) {stop("Can't estimate nugget now")}
      self$param.est <- param.est
      self$useC <- useC
      self$useGrad <- useGrad
      self$parallel <- parallel
      if (self$parallel) {
        if (parallel_cores == "detect") {
          self$parallel_cores <- parallel::detectCores()
        } else {
          self$parallel_cores <- parallel_cores
        }
      } else {self$parallel_cores <- 1}
      self$restarts <- restarts
      if (optimizer %in% c("L-BFGS-B", "BFGS", "lbfgs", "genoud")) {
        self$optimizer <- optimizer
      } else {
        stop(paste0('optimizer must be one of c("L-BFGS-B", "BFGS",',
                    ' "lbfgs, "genoud")'))
      }

      self$update_K_and_estimates() # Need to get mu_hat before starting
      # self$mu_hat <- mean(Z)
      self$fit()
      invisible(self)
    },
    # initialize_GauPr = function() {
    # },
    fit = function(X, Z) {
      self$update()
    },
    update_K_and_estimates = function () {
      # Update K, Kinv, mu_hat, and s2_hat, maybe nugget too
      self$K <- self$kernel$k(self$X) + diag(self$kernel$s2 * self$nug,
                                             self$N)
      while(T) {
        try.chol <- try(self$Kchol <- chol(self$K), silent = T)
        if (!inherits(try.chol, "try-error")) {break}
        warning("Can't Cholesky, increasing nugget #7819553")
        oldnug <- self$nug
        self$nug <- max(1e-8, 2 * self$nug)
        self$K <- self$K + diag(self$kernel$s2 * (self$nug - oldnug),
                                self$N)
        print(c(oldnug, self$nug))
      }
      self$Kinv <- chol2inv(self$Kchol)
      # self$mu_hat <- sum(self$Kinv %*% self$Z) / sum(self$Kinv)
      self$mu_hatX <- self$trend$Z(X=self$X)
      self$Kinv_Z_minus_mu_hatX <- c(self$Kinv %*% (self$Z - self$mu_hatX))
      # self$s2_hat <- c(t(self$Z - self$mu_hat) %*% self$Kinv %*%
      #                               (self$Z - self$mu_hat) / self$N)
      self$s2_hat <- self$kernel$s2
    },
    predict = function(XX, se.fit=F, covmat=F, split_speed=F) {
      self$pred(XX=XX, se.fit=se.fit, covmat=covmat,
                split_speed=split_speed)
    },
    pred = function(XX, se.fit=F, covmat=F, split_speed=F) {
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
          #se <- rep(0, length(mn)) # NEG VARS will be 0 for se,
          #  NOT SURE I WANT THIS
        }

        ni <- 40 # batch size
        Nni <- ceiling(N/ni)-1
        for (j in 0:Nni) {
          XXj <- XX[(j*ni+1):(min((j+1)*ni,N)), , drop=FALSE]
          # kxxj <- self$corr_func(XXj)
          # kx.xxj <- self$corr_func(self$X, XXj)
          predj <- self$pred_one_matrix(XX=XXj, se.fit=se.fit,
                                        covmat=covmat)
          #mn[(j*ni+1):(min((j+1)*ni,N))] <- pred_meanC(XXj, kx.xxj,
          #                                self$mu_hat, self$Kinv, self$Z)
          if (!se.fit) { # if no se.fit, just set vector
            mn[(j*ni+1):(min((j+1)*ni,N))] <- predj
          } else { # otherwise set all three from data.frame
            mn[(j*ni+1):(min((j+1)*ni,N))] <- predj$mean
            #s2j <- pred_var(XXj, kxxj, kx.xxj, self$s2_hat, self$Kinv,
            #         self$Z)
            #s2[(j*ni+1):(min((j+1)*ni,N))] <- s2j
            s2[(j*ni+1):(min((j+1)*ni,N))] <- predj$s2
            se[(j*ni+1):(min((j+1)*ni,N))] <- predj$se

          }
        }
        #se[s2>=0] <- sqrt(s2[s2>=0])

        # # Unnormalize if needed
        # if (self$normalize) {
        #   mn <- mn * self$normalize_sd + self$normalize_mean
        #   if (se.fit) {
        #     se <- se * self$normalize_sd
        #     s2 <- s2 * self$normalize_sd^2
        #   }
        # }

        if (!se.fit) {# covmat is always FALSE for split_speed } &
          # !covmat) {
          return(mn)
        } else {
          return(data.frame(mean=mn, s2=s2, se=se))
        }

      } else {
        pred1 <- self$pred_one_matrix(XX=XX, se.fit=se.fit,
                                      covmat=covmat, return_df=TRUE)
        return(pred1)
      }
    },
    pred_one_matrix = function(XX, se.fit=F, covmat=F, return_df=FALSE) {
      # input should already be checked for matrix
      # kxx <- self$kernel$k(XX) + diag(self$nug * self$s2_hat, nrow(XX))
      kx.xx <- self$kernel$k(self$X, XX)
      # mn <- pred_meanC(XX, kx.xx, self$mu_hat, self$Kinv, self$Z)
      # Changing to use trend, mu_hat is matrix
      # mu_hat_matX <- self$trend$Z(self$X)
      mu_hat_matXX <- self$trend$Z(XX)

      # mn <- pred_meanC_mumat(XX, kx.xx, self$mu_hatX, mu_hat_matXX,
      #                        self$Kinv, self$Z)
      # New way using _fast is O(n^2)
      mn <- pred_meanC_mumat_fast(XX, kx.xx, self$Kinv_Z_minus_mu_hatX,
                                  mu_hat_matXX)

      if (self$normalize) {
        mn <- mn * self$normalize_sd + self$normalize_mean
      }

      if (!se.fit & !covmat) {
        return(mn)
      }
      if (covmat) {
        # new for kernel
        kxx <- self$kernel$k(XX) + diag(self$nug * self$s2_hat, nrow(XX))
        covmatdat <- kxx - t(kx.xx) %*% self$Kinv %*% kx.xx

        if (self$normalize) {
          covmatdat <- covmatdat * self$normalize_sd ^ 2
        }

        # #covmatdat <- self$pred_var(XX, kxx=kxx, kx.xx=kx.xx, covmat=T)
        # covmatdat <- pred_cov(XX, kxx, kx.xx, self$s2_hat, self$Kinv,
        #                       self$Z)
        s2 <- diag(covmatdat)
        se <- rep(1e-8, length(mn)) # NEG VARS will be 0 for se,
        #  NOT SURE I WANT THIS
        se[s2>=0] <- sqrt(s2[s2>=0])
        return(list(mean=mn, s2=s2, se=se, cov=covmatdat))
      }


      # new for kernel
      # covmatdat <- kxx - t(kx.xx) %*% self$Kinv %*% kx.xx
      # s2 <- diag(covmatdat)
      # Better way doesn't do full matmul twice, 2x speed for 50 rows,
      #                  20x speedup for 1000 rows
      # This method is bad since only diag of k(XX) is needed
      # kxx <- self$kernel$k(XX) + diag(self$nug * self$s2_hat, nrow(XX))
      # s2 <- diag(kxx) - colSums( (kx.xx) * (self$Kinv %*% kx.xx))
      # This is bad since apply is actually really slow for a
      #                         simple function like this
      # diag.kxx <- self$nug * self$s2_hat + apply(XX, 1,
      #                 function(xrow) {self$kernel$k(xrow)})
      # s2 <- diag.kxx - colSums( (kx.xx) * (self$Kinv %*% kx.xx))
      # This method is fastest, assumes that correlation of point
      #   with itself is 1, which is true for basic kernels.
      diag.kxx <- self$nug * self$s2_hat + rep(self$s2_hat, nrow(XX))
      s2 <- diag.kxx - colSums( (kx.xx) * (self$Kinv %*% kx.xx))

      if (self$normalize) {
        s2 <- s2 * self$normalize_sd ^ 2
      }

      # s2 <- pred_var(XX, kxx, kx.xx, self$s2_hat, self$Kinv, self$Z)
      se <- rep(0, length(mn)) # NEG VARS will be 0 for se,
      #   NOT SURE I WANT THIS
      se[s2>=0] <- sqrt(s2[s2>=0])

      # se.fit but not covmat
      if (return_df) {
        # data.frame is really slow compared to cbind or list
        data.frame(mean=mn, s2=s2, se=se)
      } else {
        list(mean=mn, s2=s2, se=se)
      }
    },
    pred_mean = function(XX, kx.xx) { # 2-8x faster to use pred_meanC
      # c(self$mu_hat + t(kx.xx) %*% self$Kinv %*% (self$Z - self$mu_hat))
      # mu_hat_matX <- self$trend$Z(self$X)
      mu_hat_matXX <- self$trend$Z(XX)
      c(mu_hat_matXX + t(kx.xx) %*% self$Kinv %*% (self$Z - self$mu_hatX))
    },
    pred_meanC = function(XX, kx.xx) { # Don't use if R uses pass by copy(?)
      # pred_meanC(XX, kx.xx, self$mu_hat, self$Kinv, self$Z)
      # mu_hat_matX <- self$trend$Z(self$X)
      mu_hat_matXX <- self$trend$Z(XX)
      # This way is O(n^2)
      # pred_meanC_mumat(XX, kx.xx, self$mu_hatX, mu_hat_matXX,
      # self$Kinv, self$Z)
      # New way is O(n), but not faster in R
      # mu_hat_matXX +
      #   colSums(sweep(kx.xx, 1, self$Kinv_Z_minus_mu_hatX, `*`))
      # Rcpp code is slightly fast for small n, 2x for bigger n,
      pred_meanC_mumat_fast(XX, kx.xx, self$Kinv_Z_minus_mu_hatX,
                            mu_hat_matXX)
    },
    pred_var = function(XX, kxx, kx.xx, covmat=F) {
      # 2-4x faster to use C functions pred_var and pred_cov
      self$s2_hat * diag(kxx - t(kx.xx) %*% self$Kinv %*% kx.xx)
    },
    loglikelihood = function(...) {#mu=self$mu_hatX, s2=self$s2_hat) {
      # -.5 * (self$N*log(s2) + log(det(self$K)) +
      #          t(self$Z - mu)%*%self$Kinv%*%(self$Z - mu)/s2)
      K <- self$kernel$k(self$X, ...)
      browser()
      K_plus_Ainv_Sigbar <- K + self$Ainv_Sigbar
      -self$N/2*log(2*pi)-.5*log(det(K_plus_Ainv_Sigbar)) -
        .5*sum(self$Z * solve(K_plus_Ainv_Sigbar, self$Z))
    },
    optim = function (restarts = 5, param_update = T,
                      nug.update = self$nug.est, parallel=self$parallel,
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
      optim_functions <- self$get_optim_functions(
        param_update=param_update,
        nug.update=nug.update)
      #optim.func <- self$get_optim_func(param_update=param_update,
      #              nug.update=nug.update)
      # optim.grad <- self$get_optim_grad(param_update=param_update,
      #                                   nug.update=nug.update)
      # optim.fngr <- self$get_optim_fngr(param_update=param_update,
      #                                   nug.update=nug.update)
      optim.func <- optim_functions[[1]]
      optim.grad <- optim_functions[[2]]
      optim.fngr <- optim_functions[[3]]


      # # Set starting parameters and bounds
      # lower <- c()
      # upper <- c()
      # start.par <- c()
      # start.par0 <- c() # Some default params
      # if (param_update) {
      #   lower <- c(lower, self$param_optim_lower())
      #              #rep(-5, self$theta_length))
      #   upper <- c(upper, self$param_optim_upper())
      #                #rep(7, self$theta_length))
      #   start.par <- c(start.par, self$param_optim_start())
      #               #log(self$theta_short, 10))
      #   start.par0 <- c(start.par0, self$param_optim_start0())
      #                 #rep(0, self$theta_length))
      # }
      # if (nug.update) {
      #   lower <- c(lower, log(self$nug.min,10))
      #   upper <- c(upper, Inf)
      #   start.par <- c(start.par, log(self$nug,10))
      #   start.par0 <- c(start.par0, -4)
      # }
      #

      # Changing so all are gotten by self function
      lower <- self$param_optim_lower(nug.update=nug.update)
      upper <- self$param_optim_upper(nug.update=nug.update)
      # start.par <- self$param_optim_start(nug.update=nug.update)
      # start.par0 <- self$param_optim_start0(nug.update=nug.update)
      #
      param_optim_start_mat <- self$param_optim_start_mat(restarts=restarts,
                                                          nug.update=nug.update,
                                                          l=length(lower))
      if (!is.matrix(param_optim_start_mat)) {
        # Is a vector, should be a matrix with one row since it applies
        #    over columns
        param_optim_start_mat <- matrix(param_optim_start_mat, nrow=1)
      }

      # This will make sure it at least can start
      # Run before it sets initial parameters
      # try.devlog <- try(devlog <- optim.func(start.par), silent = T)
      try.devlog <- try(devlog <- optim.func(param_optim_start_mat[,1]),
                        silent = T)
      if (inherits(try.devlog, "try-error")) {
        warning("Current nugget doesn't work, increasing it #31973")
        # This will increase the nugget until cholesky works
        self$update_K_and_estimates()
        # devlog <- optim.func(start.par)
        devlog <- optim.func(param_optim_start_mat[,1])
      }

      # Find best params with optimization, start with current params in
      #       case all give error
      # Current params
      #best <- list(par=c(log(self$theta_short, 10), log(self$nug,10)),
      #                     value = devlog)
      # best <- list(par=start.par, value = devlog)
      best <- list(par=param_optim_start_mat[,1], value = devlog)
      if (self$verbose >= 2) {
        cat("Optimizing\n");cat("\tInitial values:\n");print(best)
      }
      #details <- data.frame(start=paste(c(self$theta_short,self$nug),
      #   collapse=","),end=NA,value=best$value,func_evals=1,
      #   grad_evals=NA,convergence=NA, message=NA, stringsAsFactors=F)
      details <- data.frame(
        start=paste(param_optim_start_mat[,1],collapse=","),end=NA,
        value=best$value,func_evals=1,grad_evals=NA,convergence=NA,
        message=NA, stringsAsFactors=F
      )


      # runs them in parallel, first starts from current,
      #         rest are jittered or random
      sys_name <- Sys.info()["sysname"]
      if (!self$parallel) {
        # Not parallel, just use lapply
        restarts.out <- lapply(
          1:(1+restarts),
          function(i){
            self$optimRestart(start.par=start.par,
                              start.par0=start.par0,
                              param_update=param_update,
                              nug.update=nug.update,
                              optim.func=optim.func,
                              optim.grad=optim.grad,
                              optim.fngr=optim.fngr,
                              lower=lower, upper=upper,
                              jit=(i!=1),
                              start.par.i=param_optim_start_mat[,i])})
      } else if (sys_name == "Windows") {
        # Parallel on Windows
        #  Not much speedup since it has to copy each time.
        #  Only maybe worth it on big problems.
        parallel_cluster <- parallel::makeCluster(
          spec = self$parallel_cores, type = "SOCK")
        restarts.out <- parallel::clusterApplyLB(
          cl=parallel_cluster,
          1:(1+restarts),
          function(i){
            self$optimRestart(start.par=start.par,
                              start.par0=start.par0,
                              param_update=param_update,
                              nug.update=nug.update,
                              optim.func=optim.func,
                              optim.grad=optim.grad,
                              optim.fngr=optim.fngr,
                              lower=lower, upper=upper,
                              jit=(i!=1),
                              start.par.i=param_optim_start_mat[,i])})
        parallel::stopCluster(parallel_cluster)
        #, mc.cores = parallel_cores)
      } else { # Mac/Unix
        restarts.out <- parallel::mclapply(
          1:(1+restarts),
          function(i){
            self$optimRestart(start.par=start.par,
                              start.par0=start.par0,
                              param_update=param_update,
                              nug.update=nug.update,
                              optim.func=optim.func,
                              optim.grad=optim.grad,
                              optim.fngr=optim.fngr,
                              lower=lower,
                              upper=upper,
                              jit=(i!=1))},
          start.par.i=param_optim_start_mat[,i],
          mc.cores = parallel_cores)
      }
      new.details <- t(sapply(restarts.out,function(dd){dd$deta}))
      vals <- sapply(restarts.out,
                     function(ii){
                       if (inherits(ii$current,"try-error")){Inf}
                       else ii$current$val
                     }
      )
      bestparallel <- which.min(vals) #which.min(new.details$value)
      if(inherits(
        try(restarts.out[[bestparallel]]$current$val, silent = T),
        "try-error")
      ) { # need this in case all are restart vals are Inf
        print("All restarts had error, keeping initial")
      } else if (restarts.out[[bestparallel]]$current$val < best$val) {
        best <- restarts.out[[bestparallel]]$current
      }
      details <- rbind(details, new.details)

      if (self$verbose >= 2) {print(details)}

      # If new nug is below nug.min, optimize again with fixed nug
      # Moved into update_params, since I don't want to set nugget here

      if (nug.update) {
        best$par[length(best$par)] <- 10 ^ (best$par[length(best$par)])
      }
      best
    },
    optimRestart = function (start.par, start.par0, param_update,
                             nug.update, optim.func, optim.grad,
                             optim.fngr, lower, upper, jit=T,
                             start.par.i) {
      #
      # FOR lognug RIGHT NOW, seems to be at least as fast,
      #    up to 5x on big data, many fewer func_evals
      #    still want to check if it is better or not

      # if (runif(1) < .33 & jit) {
      #   # restart near some spot to avoid getting stuck in bad spot
      #   start.par.i <- start.par0
      #   #print("start at zero par")
      # } else { # jitter from current params
      #   start.par.i <- start.par
      # }
      # if (FALSE) {#jit) {
      #   #if (param_update) {start.par.i[1:self$theta_length] <-
      #        start.par.i[1:self$theta_length] +
      #        rnorm(self$theta_length,0,2)} # jitter betas
      #   theta_indices <- 1:length(self$param_optim_start())
      #                          #if () -length(start.par.i)
      #   if (param_update) {start.par.i[theta_indices] <-
      #        start.par.i[theta_indices] +
      #        self$param_optim_jitter(start.par.i[theta_indices])}
      #        # jitter betas
      #   if (nug.update) {start.par.i[length(start.par.i)] <-
      #        start.par.i[length(start.par.i)] + min(4, rexp(1,1))}
      #        # jitter nugget
      # }

      # if (runif(1) < .33) { # Start at 0 params
      #   start.par.i <- self$kernel$param_optim_start0(jitter=jit)
      # } else { # Start at current params
      #   start.par.i <- self$kernel$param_optim_start(jitter=jit)
      # }
      #
      if (self$verbose >= 2) {
        cat("\tRestart (parallel): starts pars =",start.par.i,"\n")
      }
      current <- try(
        if (self$useGrad) {
          if (is.null(optim.fngr)) {
            lbfgs::lbfgs(optim.func, optim.grad, start.par.i, invisible=1)
          } else {
            # Two options for shared grad
            if (self$optimizer == "L-BFGS-B") {
              # optim uses L-BFGS-B which uses upper and lower
              optim_share(fngr=optim.fngr, par=start.par.i,
                          method='L-BFGS-B', upper=upper, lower=lower)
            } else if (self$optimizer == "lbfgs") {
              # lbfgs does not, so no longer using it
              lbfgs_share(optim.fngr, start.par.i, invisible=1)
              # 1.7x speedup uses grad_share
            } else if (self$optimizer == "genoud") {
              capture.output(suppressWarnings({
                tmp <- rgenoud::genoud(fn=optim.func,
                                       nvars=length(start.par.i),
                                       starting.values=start.par.i,
                                       Domains=cbind(lower, upper),
                                       gr=optim.grad,
                                       boundary.enforcement = 2,
                                       pop.size=1e2, max.generations=10)
              }))
              tmp
            } else{
              stop("Optimizer not recognized")
            }
          }
        } else {
          optim(start.par.i, optim.func, method="L-BFGS-B",
                lower=lower, upper=upper, hessian=F)
        }
      )
      if (!inherits(current, "try-error")) {
        if (self$useGrad) {
          current$counts <- c(NA,NA)
          if(is.null(current$message)) {current$message=NA}
        }
        details.new <- data.frame(
          start=paste(signif(start.par.i,3),collapse=","),
          end=paste(signif(current$par,3),collapse=","),
          value=current$value,func_evals=current$counts[1],
          grad_evals=current$counts[2],
          convergence=if (is.null(current$convergence)) {NA}
          else {current$convergence},
          message=current$message, row.names = NULL, stringsAsFactors=F
        )
      } else{
        details.new <- data.frame(
          start=paste(signif(start.par.i,3),collapse=","),
          end="try-error",value=NA,func_evals=NA,grad_evals=NA,
          convergence=NA, message=current[1], stringsAsFactors=F)
      }
      list(current=current, details=details.new)
    },
    update = function (Xnew=NULL, Znew=NULL, Xall=NULL, Zall=NULL,
                       restarts = self$restarts,
                       param_update = self$param.est,
                       nug.update = self$nug.est, no_update=FALSE) {
      # Doesn't update Kinv, etc
      self$update_data(Xnew=Xnew, Znew=Znew, Xall=Xall, Zall=Zall)

      if (!no_update && (param_update || nug.update)) {
        # This option lets it skip parameter optimization entirely
        self$update_params(restarts=restarts,
                           param_update=param_update,
                           nug.update=nug.update)
      }

      self$update_K_and_estimates()

      invisible(self)
    }
  )
)
