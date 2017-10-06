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
#' @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' n <- 12
#' x <- matrix(seq(0,1,length.out = n), ncol=1)
#' y <- sin(2*pi*x) + rnorm(n,0,1e-1)
#' gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian$new(1), parallel=FALSE)
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
#' @field parallel_cores How many cores are there? It will self detect, do not set yourself.
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/lightning-viz/lightining-r/}
#'   \item{\code{new(X, Z, corr="Gauss", verbose=0, separable=T, useC=F,useGrad=T,
#'          parallel=T, nug.est=T, ...)}}{This method is used to create object of this class with \code{X} and \code{Z} as the data.}
#'
#'   \item{\code{update(Xnew=NULL, Znew=NULL, Xall=NULL, Zall=NULL,
#' restarts = 5,
#' param_update = T, nug.update = self$nug.est)}}{This method updates the model, adding new data if given, then running optimization again.}
#'   }
GauPro_kernel_model <- R6::R6Class(classname = "GauPro",
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
        param.est = NULL, # Whether parameters besides nugget (theta) should be updated
        # mu_hat = NULL,
        mu_hatX = NULL,
        s2_hat = NULL,
        # corr_func = function(...){}, # When this was NULL the child didn't overwrite with own method, it stayed as NULL
        K = NULL,
        Kchol = NULL,
        Kinv = NULL,
        verbose = 0,
        useC = TRUE,
        useGrad = FALSE,
        parallel = NULL,
        parallel_cores = NULL,
        restarts = NULL,
        normalize = NULL, # Should the Z values be normalized for internal computations?
        normalize_mean = NULL,
        normalize_sd = NULL,
        #deviance_out = NULL, #(theta, nug)
        #deviance_grad_out = NULL, #(theta, nug, overwhat)
        #deviance_fngr_out = NULL,
        initialize = function(X, Z,
                              kernel, trend,
                              verbose=0, useC=F,useGrad=T,
                              parallel=FALSE,
                              nug=1e-6, nug.min=1e-8, nug.max=Inf, nug.est=TRUE,
                              param.est = TRUE, restarts = 5,
                              normalize = FALSE,
                              ...) {
          #self$initialize_GauPr(X=X,Z=Z,verbose=verbose,useC=useC,useGrad=useGrad,
          #                      parallel=parallel, nug.est=nug.est)
          self$X <- X
          self$Z <- matrix(Z, ncol=1)
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
          if ("R6ClassGenerator" %in% class(kernel)) { # Let generator be given so D can be set auto
            self$kernel <- kernel$new(D=self$D)
          } else if ("GauPro_kernel" %in% class(kernel)) { # Otherwise it should already be a kernel
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

          self$nug <- nug
          self$nug.min <- nug.min
          self$nug.max <- nug.max
          self$nug.est <- nug.est
          # if (nug.est) {stop("Can't estimate nugget now")}
          self$param.est <- param.est
          self$useC <- useC
          self$useGrad <- useGrad
          self$parallel <- parallel
          if (self$parallel) {self$parallel_cores <- parallel::detectCores()}
          else {self$parallel_cores <- 1}
          self$restarts <- restarts

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
          self$K <- self$kernel$k(self$X) + diag(self$kernel$s2 * self$nug, self$N)
          while(T) {
            try.chol <- try(self$Kchol <- chol(self$K), silent = T)
            if (!inherits(try.chol, "try-error")) {break}
            warning("Can't Cholesky, increasing nugget #7819553")
            oldnug <- self$nug
            self$nug <- max(1e-8, 2 * self$nug)
            self$K <- self$K + diag(self$kernel$s2 * (self$nug - oldnug), self$N)
            print(c(oldnug, self$nug))
          }
          self$Kinv <- chol2inv(self$Kchol)
          # self$mu_hat <- sum(self$Kinv %*% self$Z) / sum(self$Kinv)
          self$mu_hatX <- self$trend$Z(X=self$X)
          # self$s2_hat <- c(t(self$Z - self$mu_hat) %*% self$Kinv %*% (self$Z - self$mu_hat) / self$N)
          self$s2_hat <- self$kernel$s2
        },
        predict = function(XX, se.fit=F, covmat=F, split_speed=T) {
          self$pred(XX=XX, se.fit=se.fit, covmat=covmat, split_speed=split_speed)
        },
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

            # # Unnormalize if needed
            # if (self$normalize) {
            #   mn <- mn * self$normalize_sd + self$normalize_mean
            #   if (se.fit) {
            #     se <- se * self$normalize_sd
            #     s2 <- s2 * self$normalize_sd^2
            #   }
            # }

            if (!se.fit) {# covmat is always FALSE for split_speed } & !covmat) {
              return(mn)
            } else {
              return(data.frame(mean=mn, s2=s2, se=se))
            }

          } else {
            pred1 <- self$pred_one_matrix(XX=XX, se.fit=se.fit, covmat=covmat)
            return(pred1)
          }
        },
        pred_one_matrix = function(XX, se.fit=F, covmat=F) {
          # input should already be check for matrix
          kxx <- self$kernel$k(XX) + self$nug
          kx.xx <- self$kernel$k(self$X, XX)
          # mn <- pred_meanC(XX, kx.xx, self$mu_hat, self$Kinv, self$Z)
          # Changing to use trend, mu_hat is matrix
          # mu_hat_matX <- self$trend$Z(self$X)
          mu_hat_matXX <- self$trend$Z(XX)

          mn <- pred_meanC_mumat(XX, kx.xx, self$mu_hatX, mu_hat_matXX, self$Kinv, self$Z)

          if (self$normalize) {mn <- mn * self$normalize_sd + self$normalize_mean}

          if (!se.fit & !covmat) {
            return(mn)
          }
          if (covmat) {
            # new for kernel
            covmatdat <- kxx - t(kx.xx) %*% self$Kinv %*% kx.xx

            if (self$normalize) {
              covmatdat <- covmatdat * self$normalize_sd ^ 2
            }

            # #covmatdat <- self$pred_var(XX, kxx=kxx, kx.xx=kx.xx, covmat=T)
            # covmatdat <- pred_cov(XX, kxx, kx.xx, self$s2_hat, self$Kinv, self$Z)
            s2 <- diag(covmatdat)
            se <- rep(1e-8, length(mn)) # NEG VARS will be 0 for se, NOT SURE I WANT THIS
            se[s2>=0] <- sqrt(s2[s2>=0])
            return(list(mean=mn, s2=s2, se=se, cov=covmatdat))
          }


          # new for kernel
          covmatdat <- kxx - t(kx.xx) %*% self$Kinv %*% kx.xx
          s2 <- diag(covmatdat)

          if (self$normalize) {
            s2 <- s2 * self$normalize_sd ^ 2
          }

          # s2 <- pred_var(XX, kxx, kx.xx, self$s2_hat, self$Kinv, self$Z)
          se <- rep(0, length(mn)) # NEG VARS will be 0 for se, NOT SURE I WANT THIS
          se[s2>=0] <- sqrt(s2[s2>=0])

          # se.fit but not covmat
          data.frame(mean=mn, s2=s2, se=se)
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
          pred_meanC_mumat(XX, kx.xx, self$mu_hatX, mu_hat_matXX, self$Kinv, self$Z)
        },
        pred_var = function(XX, kxx, kx.xx, covmat=F) { # 2-4x faster to use C functions pred_var and pred_cov
          self$s2_hat * diag(kxx - t(kx.xx) %*% self$Kinv %*% kx.xx)
        },
        pred_LOO = function(se.fit=FALSE) {#browser()
          # Predict LOO (leave-one-out) on data used to fit model
          # See vignette for explanation of equations
          # If se.fit==T, then calculate the LOO se and the corresponding t score
          Z_LOO <- numeric(self$N)
          if (se.fit) {Z_LOO_s2 <- numeric(self$N)}
          Z_trend <- self$trend$Z(self$X)
          for (i in 1:self$N) {
            E <- self$Kinv[-i, -i] # Kinv without i
            b <- self$K[    i, -i] # K    between i and rest
            g <- self$Kinv[ i, -i] # Kinv between i and rest
            Ainv <- E + E %*% b %*% g / (1-sum(g*b)) # Kinv for K if i wasn't in K
            Zi_LOO <- Z_trend[i] + c(b %*% Ainv %*% (self$Z[-i] - Z_trend[-i]))
            Z_LOO[i] <- Zi_LOO
            if (se.fit) {
              Zi_LOO_s2 <- self$K[i,i] - c(b %*% Ainv %*% b)
              Z_LOO_s2[i] <- Zi_LOO_s2
            }
          }
          if (se.fit) { # Return df with se and t if se.fit
            t_LOO <- (self$Z - Z_LOO) / Z_LOO_s2
            data.frame(fit=Z_LOO, se.fit=Z_LOO_s2, t=t_LOO)
          } else { # Else just mean LOO
            Z_LOO
          }
        },
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
          Sigma.try <- try(newy <- MASS::mvrnorm(n=n2, mu=px$mean, Sigma=px$cov), silent = TRUE)
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
          points(self$X,
                 if (self$normalize) {self$Z * self$normalize_sd + self$normalize_mean}
                   else {self$Z},
                 pch=19, col=1, cex=2)
        },
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
                 if (self$normalize) {self$Z * self$normalize_sd + self$normalize_mean}
                 else {self$Z},
                 pch=19, col=1, cex=2)
        },
        plot2D = function() {
          if (self$D != 2) {stop("plot2D only works in 2D")}
          mins <- apply(self$X, 2, min)
          maxs <- apply(self$X, 2, max)
          xmin <- mins[1] - .03 * (maxs[1] - mins[1])
          xmax <- maxs[1] + .03 * (maxs[1] - mins[1])
          ymin <- mins[2] - .03 * (maxs[2] - mins[2])
          ymax <- maxs[2] + .03 * (maxs[2] - mins[2])
          ContourFunctions::cf_func(self$predict, batchmax=Inf,
                                    xlim=c(xmin, xmax),
                                    ylim=c(ymin, ymax),
                                    pts=self$X)
        },
        loglikelihood = function(mu=self$mu_hatX, s2=self$s2_hat) {
          -.5 * (self$N*log(s2) + log(det(self$K)) + t(self$Z - mu)%*%self$Kinv%*%(self$Z - mu)/s2)
        },
        get_optim_functions = function(param_update, nug.update) {
          # self$kernel$get_optim_functions(param_update=param_update)
          # if (nug.update) { # Nug will be last in vector of parameters
          #   list(
          #     fn=function(params) {
          #       l <- length(params)
          #       self$deviance(params=params[1:(l-1)], nuglog=params[l])
          #     },
          #     gr=function(params) {
          #       l <- length(params)
          #       self$deviance_grad(params=params[1:(l-1)], nuglog=params[l], nug.update=nug.update)
          #     },
          #     fngr=function(params) {
          #       l <- length(params)
          #       list(
          #         fn=function(params) {
          #           self$deviance(params=params[1:(l-1)], nuglog=params[l])
          #         },
          #         gr=function(params) {self$deviance_grad(params=params[1:(l-1)], nuglog=params[l], nug.update=nug.update)
          #         }
          #       )
          #     }
          #   )
          # } else {
          #   list(
          #     fn=function(params) {self$deviance(params=params)},
          #     gr=function(params) {self$deviance_grad(params=params, nug.update=nug.update)},
          #     fngr=function(params) {
          #       list(
          #         fn=function(params) {self$deviance(params=params)},
          #         gr=function(params) {self$deviance_grad(params=params, nug.update=nug.update)}
          #       )
          #     }
          #   )
          # }


          # Need to get trend, kernel, and nug separated out into args
          tl <- length(self$trend$param_optim_start())
          kl <- length(self$kernel$param_optim_start())
          nl <- as.integer(nug.update)
          ti <- if (tl>0) {1:tl} else {c()}
          ki <- if (kl>0) {tl + 1:kl} else {c()}
          ni <- if (nl>0) {tl+kl+1:nl} else {c()}
          list(
            fn=function(params) {
              tparams <- if (tl>0) {params[ti]} else {NULL}
              kparams <- if (kl>0) {params[ki]} else {NULL}
              nparams <- if (nl>0) {params[ni]} else {NULL}
              self$deviance(params=kparams, nuglog=nparams, trend_params=tparams)
            },
            gr=function(params) {
              tparams <- if (tl>0) {params[ti]} else {NULL}
              kparams <- if (kl>0) {params[ki]} else {NULL}
              nparams <- if (nl>0) {params[ni]} else {NULL}
              self$deviance_grad(params=kparams, nuglog=nparams, trend_params=tparams, nug.update=nug.update)
            },
            fngr=function(params) {
              tparams <- if (tl>0) {params[ti]} else {NULL}
              kparams <- if (kl>0) {params[ki]} else {NULL}
              nparams <- if (nl>0) {params[ni]} else {NULL}
              self$deviance_fngr(params=kparams, nuglog=nparams, trend_params=tparams, nug.update=nug.update)
            }
          )
        },
        param_optim_lower = function(nug.update) {
          if (nug.update) {
          #   c(self$kernel$param_optim_lower(), log(self$nug.min,10))
            nug_lower <- log(self$nug.min, 10)
          } else {
          #   self$kernel$param_optim_lower()
            nug_lower <- c()
          }
          trend_lower <- self$trend$param_optim_lower()
          kern_lower <- self$kernel$param_optim_lower()
          c(trend_lower, kern_lower, nug_lower)
        },
        param_optim_upper = function(nug.update) {
          if (nug.update) {
          #   c(self$kernel$param_optim_upper(), Inf)
            nug_upper <- log(self$nug.max, 10)
          } else {
          #   self$kernel$param_optim_upper()
            nug_upper <- c()
          }
          trend_upper <- self$trend$param_optim_upper()
          kern_upper <- self$kernel$param_optim_upper()
          c(trend_upper, kern_upper, nug_upper)

        },
        param_optim_start = function(nug.update, jitter) {
          # param_start <- self$kernel$param_optim_start(jitter=jitter)
          if (nug.update) {
            nug_start <- log(self$nug,10)
            if (jitter) {nug_start <- nug_start + rexp(1, 1)}
            # c(param_start, nug_start)
          } else {
            # param_start
            nug_start <- c()
          }

          trend_start <- self$trend$param_optim_start()
          kern_start <- self$kernel$param_optim_start(jitter=jitter)
          # nug_start <- Inf
          c(trend_start, kern_start, nug_start)
        },
        param_optim_start0 = function(nug.update, jitter) {
          # param_start <- self$kernel$param_optim_start0(jitter=jitter)
          if (nug.update) {
            nug_start <- -4
            if (jitter) {nug_start <- nug_start + rexp(1, 1)}
            # c(param_start, nug_start)
          } else {
            # param_start
            nug_start <- c()
          }
          trend_start <- self$trend$param_optim_start()
          kern_start <- self$kernel$param_optim_start(jitter=jitter)
          # nug_start <- Inf
          c(trend_start, kern_start, nug_start)
        },
        param_optim_start_mat = function(restarts, nug.update, l) {
          s0 <- sample(c(T,F), size=restarts+1, replace=TRUE, prob = c(.33,.67))
          s0[1] <- TRUE
          sapply(1:(restarts+1), function(i) {
            if (s0[i]) {
              self$param_optim_start0(nug.update=nug.update, jitter=(i!=1))
            } else {
              self$param_optim_start(nug.update=nug.update, jitter=(i!=1))
            }
          })
          # mat <- matrix(0, nrow=restarts, ncol=l)
          # mat[1,] <- self$param_optim_start0(nug.update=nug.update)
        },
        optim = function (restarts = 5, param_update = T, nug.update = self$nug.est, parallel=self$parallel, parallel_cores=self$parallel_cores) {
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


          # # Set starting parameters and bounds
          # lower <- c()
          # upper <- c()
          # start.par <- c()
          # start.par0 <- c() # Some default params
          # if (param_update) {
          #   lower <- c(lower, self$param_optim_lower())#rep(-5, self$theta_length))
          #   upper <- c(upper, self$param_optim_upper())#rep(7, self$theta_length))
          #   start.par <- c(start.par, self$param_optim_start())#log(self$theta_short, 10))
          #   start.par0 <- c(start.par0, self$param_optim_start0())#rep(0, self$theta_length))
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


          # This will make sure it at least can start
          # Run before it sets initial parameters
          # try.devlog <- try(devlog <- optim.func(start.par), silent = T)
          try.devlog <- try(devlog <- optim.func(param_optim_start_mat[,1]), silent = T)
          if (inherits(try.devlog, "try-error")) {
            warning("Current nugget doesn't work, increasing it #31973")
            self$update_K_and_estimates() # This will increase the nugget until cholesky works
            # devlog <- optim.func(start.par)
            devlog <- optim.func(param_optim_start_mat[,1])
          }

          # Find best params with optimization, start with current params in case all give error
          # Current params
          #best <- list(par=c(log(self$theta_short, 10), log(self$nug,10)), value = devlog)
          # best <- list(par=start.par, value = devlog)
          best <- list(par=param_optim_start_mat[,1], value = devlog)
          if (self$verbose >= 2) {cat("Optimizing\n");cat("\tInitial values:\n");print(best)}
          #details <- data.frame(start=paste(c(self$theta_short,self$nug),collapse=","),end=NA,value=best$value,func_evals=1,grad_evals=NA,convergence=NA, message=NA, stringsAsFactors=F)
          details <- data.frame(start=paste(param_optim_start_mat[,1],collapse=","),end=NA,value=best$value,func_evals=1,grad_evals=NA,convergence=NA, message=NA, stringsAsFactors=F)


          # runs them in parallel, first starts from current, rest are jittered or random
          sys_name <- Sys.info()["sysname"]
          if (sys_name == "Windows" | !self$parallel) {
            # Trying this so it works on Windows
            restarts.out <- lapply( 1:(1+restarts), function(i){self$optimRestart(start.par=start.par, start.par0=start.par0, param_update=param_update, nug.update=nug.update, optim.func=optim.func, optim.grad=optim.grad, optim.fngr=optim.fngr, lower=lower, upper=upper, jit=(i!=1), start.par.i=param_optim_start_mat[,i])})#, mc.cores = parallel_cores)
          } else { # Mac/Unix
            restarts.out <- parallel::mclapply(1:(1+restarts), function(i){self$optimRestart(start.par=start.par, start.par0=start.par0, param_update=param_update, nug.update=nug.update, optim.func=optim.func, optim.grad=optim.grad, optim.fngr=optim.fngr,lower=lower, upper=upper, jit=(i!=1))}, start.par.i=param_optim_start_mat[,i], mc.cores = parallel_cores)
          }
          new.details <- t(sapply(restarts.out,function(dd){dd$deta}))
          vals <- sapply(restarts.out,
                         function(ii){
                           if (inherits(ii$current,"try-error")){Inf}
                           else ii$current$val
                         }
          )
          bestparallel <- which.min(vals) #which.min(new.details$value)
          if(inherits(try(restarts.out[[bestparallel]]$current$val, silent = T), "try-error")) { # need this in case all are restart vals are Inf
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
        optimRestart = function (start.par, start.par0, param_update, nug.update, optim.func, optim.grad, optim.fngr, lower, upper, jit=T, start.par.i) {
          #
          # FOR lognug RIGHT NOW, seems to be at least as fast, up to 5x on big data, many fewer func_evals
          #    still want to check if it is better or not

          # if (runif(1) < .33 & jit) { # restart near some spot to avoid getting stuck in bad spot
          #   start.par.i <- start.par0
          #   #print("start at zero par")
          # } else { # jitter from current params
          #   start.par.i <- start.par
          # }
          # if (FALSE) {#jit) {
          #   #if (param_update) {start.par.i[1:self$theta_length] <- start.par.i[1:self$theta_length] + rnorm(self$theta_length,0,2)} # jitter betas
          #   theta_indices <- 1:length(self$param_optim_start()) #if () -length(start.par.i)
          #   if (param_update) {start.par.i[theta_indices] <- start.par.i[theta_indices] + self$param_optim_jitter(start.par.i[theta_indices])} # jitter betas
          #   if (nug.update) {start.par.i[length(start.par.i)] <- start.par.i[length(start.par.i)] + min(4, rexp(1,1))} # jitter nugget
          # }

          # if (runif(1) < .33) { # Start at 0 params
          #   start.par.i <- self$kernel$param_optim_start0(jitter=jit)
          # } else { # Start at current params
          #   start.par.i <- self$kernel$param_optim_start(jitter=jit)
          # }
          #
          if (self$verbose >= 2) {cat("\tRestart (parallel): starts pars =",start.par.i,"\n")}
          current <- try(
            if (self$useGrad) {
              if (is.null(optim.fngr)) {
                lbfgs::lbfgs(optim.func, optim.grad, start.par.i, invisible=1)
              } else {
                # Two options for shared grad
                if (TRUE) { # optim uses L-BFGS-B which uses upper and lower
                  optim_share(fngr=optim.fngr, par=start.par.i, method='L-BFGS-B', upper=upper, lower=lower)
                } else { # lbfgs does not, so no longer using it
                  lbfgs_share(optim.fngr, start.par.i, invisible=1) # 1.7x speedup uses grad_share
                }
              }
            } else {
              optim(start.par.i, optim.func, method="L-BFGS-B", lower=lower, upper=upper, hessian=F)
            }
          )
          if (!inherits(current, "try-error")) {
            if (self$useGrad) {current$counts <- c(NA,NA);if(is.null(current$message))current$message=NA}
            details.new <- data.frame(start=paste(signif(start.par.i,3),collapse=","),end=paste(signif(current$par,3),collapse=","),value=current$value,func_evals=current$counts[1],grad_evals=current$counts[2],convergence=current$convergence, message=current$message, row.names = NULL, stringsAsFactors=F)
          } else{
            details.new <- data.frame(start=paste(signif(start.par.i,3),collapse=","),end="try-error",value=NA,func_evals=NA,grad_evals=NA,convergence=NA, message=current[1], stringsAsFactors=F)
          }
          list(current=current, details=details.new)
        },
        update = function (Xnew=NULL, Znew=NULL, Xall=NULL, Zall=NULL,
                           restarts = self$restarts,
                           param_update = self$param.est, nug.update = self$nug.est, no_update=FALSE) {
          self$update_data(Xnew=Xnew, Znew=Znew, Xall=Xall, Zall=Zall) # Doesn't update Kinv, etc

          if (!no_update && (param_update || nug.update)) { # This option lets it skip parameter optimization entirely
            self$update_params(restarts=restarts, param_update=param_update,nug.update=nug.update)
          }

          self$update_K_and_estimates()

          invisible(self)
        },
        update_params = function(..., nug.update) {
          # start_params = self$kernel$get_optim_start_params()
          optim_out <- self$optim(..., nug.update=nug.update)
          # lpar <- length(optim_out$par)
          tl <- length(self$trend$param_optim_start())
          kl <- length(self$kernel$param_optim_start())
          nl <- as.integer(nug.update)
          ti <- if (tl>0) {1:tl} else {c()}
          ki <- if (kl>0) {tl + 1:kl} else {c()}
          ni <- if (nl>0) {tl+kl+1:nl} else {c()}
          if (nug.update) {
            # self$nug <- optim_out$par[lpar] # optim already does 10^
            # self$kernel$set_params_from_optim(optim_out$par[1:(lpar-1)])
            self$nug <- optim_out$par[ni] # optim already does 10^
          } else {
            # self$kernel$set_params_from_optim(optim_out$par)
          }
          self$kernel$set_params_from_optim(optim_out$par[ki])
          self$trend$set_params_from_optim(optim_out$par[ti])
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
            if (self$normalize) {
              self$normalize_mean <- mean(self$Z)
              self$normalize_sd <- sd(self$Z)
              self$Z <- (self$Z - self$normalize_mean) / self$normalize_sd
            }
          } else if (!is.null(Znew)) {
            Znewmat <- if (is.matrix(Znew)) Znew else matrix(Znew,ncol=1)
            if (self$normalize) {Znewmat <- (Znewmat - self$normalize_mean) / self$normalize_sd}
            self$Z <- rbind(self$Z, Znewmat)
          }
          #if (!is.null(Xall) | !is.null(Xnew)) {self$update_K_and_estimates()} # update Kinv, etc, DONT THINK I NEED IT
        },
        update_corrparams = function (...) {
          self$update(nug.update = F, ...=...)
        },
        update_nugget = function (...) {
          self$update(param_update = F, ...=...)
        },
        # deviance_searchnug = function() {
        #   optim(self$nug, function(nnug) {self$deviance(nug=nnug)}, method="L-BFGS-B", lower=0, upper=Inf, hessian=F)$par
        # },
        # nugget_update = function () {
        #   nug <- self$deviance_searchnug()
        #   self$nug <- nug
        #   self$update_K_and_estimates()
        # },
        deviance = function(params=NULL, nug=self$nug, nuglog, trend_params=NULL) {#browser()#print(c(params, nuglog))
          if (!missing(nuglog) && !is.null(nuglog)) {
            nug <- 10^nuglog
          }
          if (any(is.nan(params), is.nan(nug))) {if (self$verbose >= 2) {print("In deviance, returning Inf #92387")};return(Inf)}
          K <- self$kernel$k(x=self$X, params=params) +
            diag(nug, self$N) * self$kernel$s2_from_params(params=params)
          if (is.nan(log(det(K)))) {browser();return(Inf)}
          Z_hat <- self$trend$Z(X=self$X, params=trend_params)
          # dev.try <- try(dev <- log(det(K)) + sum((self$Z - self$mu_hat) * solve(K, self$Z - self$mu_hat)))
          dev.try <- try(dev <- log(det(K)) + sum((self$Z - Z_hat) * solve(K, self$Z - Z_hat)))
          if (inherits(dev.try, "try-error")) {if (self$verbose>=2) {print("Deviance error #87126, returning Inf")}; return(Inf)}
          # print(c(params, nuglog, dev))
          if (is.infinite(abs(dev))) {if (self$verbose>=2) {print("Deviance infinite #2332, returning Inf")};return(Inf)}
          dev
        },
        deviance_grad = function(params=NULL, kernel_update=TRUE,
                                 X=self$X,
                                 nug=self$nug, nug.update, nuglog,
                                 trend_params=NULL, trend_update=TRUE) {
          if (!missing(nuglog) && !is.null(nuglog)) {
            nug <- 10^nuglog
          }
          if (any(is.nan(params), is.nan(nug))) {if (self$verbose>=2) {print("In deviance_grad, returning NaN #92387")};return(rep(NaN, length(params)+as.integer(isTRUE(nug.update))))}
          C_nonug <- self$kernel$k(x=X, params=params)
          s2_from_kernel <- self$kernel$s2_from_params(params=params)
          C <- C_nonug + s2_from_kernel * diag(nug, self$N)
          dC_dparams_out <- self$kernel$dC_dparams(params=params, X=X, C=C, C_nonug=C_nonug, nug=nug)
          dC_dparams <- dC_dparams_out#[[1]] # First of list should be list of dC_dparams
          # s2_from_kernel <- dC_dparams_out[[2]] # Second should be s2 for nugget deriv
          Z_hat <- self$trend$Z(X=X, params=trend_params)
          dZ_dparams <- self$trend$dZ_dparams(X=X, params=trend_params)
          # yminusmu <- self$Z - self$mu_hat
          yminusmu <- self$Z - Z_hat
          solve.try <- try(Cinv_yminusmu <- solve(C, yminusmu))
          if (inherits(solve.try, "try-error")) { if (self$verbose>=2) {print("Deviance grad error #63466, returning Inf")};  return(Inf)}


          out <- c()
          if (length(dZ_dparams) > 0 && trend_update) {
            trend_gradfunc <- function(di) {
              -2 * t(yminusmu) %*% solve(C, di) # Siginv %*% du/db
            }
            trend_out <- apply(dZ_dparams, 2, trend_gradfunc)
            out <- trend_out
          } else {
            # trend_out <- c()
          }

          gradfunc <- function(di) {
            t1 <- sum(diag(solve(C, di)))
            t2 <- sum(Cinv_yminusmu * (di %*% Cinv_yminusmu))
            t1 - t2
          }

          if (kernel_update) {
            kernel_out <- apply(dC_dparams, 1, gradfunc)
            out <- c(out, kernel_out)
          }
          # # out <- c(sapply(dC_dparams[[1]],gradfunc), gradfunc(dC_dparams[[2]]))
          # out <- sapply(dC_dparams,gradfunc)
          if (nug.update) {
            nug_out <- gradfunc(diag(s2_from_kernel*nug*log(10), nrow(C)))
            out <- c(out, nug_out)
            # out <- c(out, gradfunc(diag(s2_from_kernel*, nrow(C)))*nug*log(10))
          }
          # print(c(params, nuglog, out))
          out
        },
        deviance_fngr = function(params=NULL, kernel_update=TRUE,
                                 X=self$X,
                                 nug=self$nug, nug.update, nuglog,
                                 trend_params=NULL, trend_update=TRUE) {#browser()
          if (!missing(nuglog) && !is.null(nuglog)) {
            nug <- 10^nuglog
          }
          if (any(is.nan(params), is.nan(nug))) {if (self$verbose>=2) {print("In deviance_grad, returning NaN #92387")};return(rep(NaN, length(params)+as.integer(isTRUE(nug.update))))}
          # C_nonug <- self$kernel$k(x=X, params=params)
          # C <- C_nonug + s2_from_kernel * diag(nug, self$N)

          # s2_from_kernel <- self$kernel$s2_from_params(params=params)
          C_dC_try <- try(
            C_dC_dparams_out <- self$kernel$C_dC_dparams(params=params, X=X, nug=nug), #C=C, C_nonug=C_nonug)
            silent = TRUE
          )
          if (inherits(C_dC_try, 'try-error')) {
            return(list(fn=self$deviance(params=params, nug=nug),
                        gr=self$deviance_grad(params=params, X=X, nug=nug, nug.update=nug.update)))
          }
          if (length(C_dC_dparams_out) < 2) {stop("Error #532987")}
          C <- C_dC_dparams_out[[1]]
          dC_dparams <- C_dC_dparams_out[[2]] # First of list should be list of dC_dparams
          # s2_from_kernel <- dC_dparams_out[[2]] # Second should be s2 for nugget deriv
          Z_hat <- self$trend$Z(X=X, params=trend_params)
          dZ_dparams <- self$trend$dZ_dparams(X=X, params=trend_params)
          # yminusmu <- self$Z - self$mu_hat
          yminusmu <- self$Z - Z_hat
          s2_from_kernel <- self$kernel$s2_from_params(params=params)
          solve.try <- try(Cinv_yminusmu <- solve(C, yminusmu))
          if (inherits(solve.try, "try-error")) { if (self$verbose>=2) {print("Deviance grad error #63466, returning Inf")};  return(Inf)}

          gr <- c()
          if (length(dZ_dparams) > 0 && trend_update) {
            trend_gradfunc <- function(di) {
              -2 * t(yminusmu) %*% solve(C, di) # Siginv %*% du/db
            }
            trend_gr <- apply(dZ_dparams, 2, trend_gradfunc)
            gr <- trend_gr
          } else {
            # trend_out <- c()
          }

          gradfunc <- function(di) {
            t1 <- sum(diag(solve(C, di)))
            t2 <- sum(Cinv_yminusmu * (di %*% Cinv_yminusmu))
            t1 - t2
          }
          # out <- c(sapply(dC_dparams[[1]],gradfunc), gradfunc(dC_dparams[[2]]))
          # gr <- sapply(dC_dparams,gradfunc)
          if (kernel_update) {
            kernel_gr <- apply(dC_dparams, 1, gradfunc)
            gr <- c(gr, kernel_gr)
          }
          if (nug.update) {
            gr <- c(gr, gradfunc(diag(s2_from_kernel*nug*log(10), nrow(C))))
            # out <- c(out, gradfunc(diag(s2_from_kernel*, nrow(C)))*nug*log(10))
          }

          # Calculate fn
          logdetC <- log(det(C))
          if (is.nan(logdetC)) {
            dev <- Inf #return(Inf)
          } else {
            dev.try <- try(dev <- logdetC + sum((yminusmu) * solve(C, yminusmu)))
            if (inherits(dev.try, "try-error")) {
              if (self$verbose>=2) {
                print("Deviance error #87126, returning Inf")
              }
              dev <- Inf #return(Inf)
            }
            # print(c(params, nuglog, dev))
            if (is.infinite(abs(dev))) {
              if (self$verbose>=2) {
                # print("Deviance infinite #2333, returning Inf")
                print("Deviance infinite #2333, returning 1e100, this is a hack and gives noticeable worse results on this restart.")
              }
              dev <- 1e100 # .Machine$double.xmax # Inf
            }
          }
          # dev

          # print(c(params, nuglog, out))
          out <- list(fn=dev, gr=gr)
          out
        },
        grad = function(XX, X=self$X, Z=self$Z) {
          if (!is.matrix(XX)) {
            if (length(XX) == self$D) {
              XX <- matrix(XX, nrow=1)
            } else {
              stop("XX must have length D or be matrix with D columns")
            }
          }
          dtrend_dx <- self$trend$dZ_dx(X=XX)
          dC_dx <- self$kernel$dC_dx(XX=XX, X=X)
          trendX <- self$trend$Z(X=X)
          Cinv_Z_minus_Zhat <- solve(self$K, Z - trendX)
          t2 <- apply(dC_dx, 1, function(U) {U %*% Cinv_Z_minus_Zhat})
          if (ncol(dtrend_dx) > 1) {
            dtrend_dx + t(t2)
          } else {
            dtrend_dx + t2
          }
          # dtrend_dx + dC_dx %*% solve(self$K, Z - trendX)
        },
        grad_norm = function (XX) {
          grad1 <- self$grad(XX)
          if (!is.matrix(grad1)) return(abs(grad1))
          apply(grad1,1, function(xx) {sqrt(sum(xx^2))})
        },
        grad_dist = function(XX) {#browser()
          nn <- nrow(XX)
          d <- ncol(XX) # or self$D
          mn <- self$grad(XX=XX)
          c2 <- self$kernel$d2C_dudv(XX=XX, X=XX)
          c1 <- self$kernel$dC_dx(XX=XX, X=self$X)
          # cv <- c2 - t(c1) %*% solve(self$Kinv, c1)
          cv <- array(data = NA, dim = c(nn, d, d))
          for (i in 1:nn) {
            tc1i <- c1[i,,] # 1D gives problem, only need transpose if D>1
            if (!is.null(dim(tc1i))) {tc1i <- t(tc1i)}
            cv[i, , ] <- c2[i,,,i] - c1[i,,] %*% (self$Kinv %*% tc1i)
          }
          list(mean=mn, cov=cv)
        },
        grad_sample = function(XX, n) {
          if (!is.matrix(XX)) {
            if (length(XX) == self$D) {XX <- matrix(XX, nrow=1)}
            else {stop("Wrong dimensions #12574")}
          }
          # if (nrow(XX) > 1) {return(apply(XX, 1, self$grad_sample))}
          if (nrow(XX) > 1) {stop("Only can do 1 grad sample at a time")}
          grad_dist <- self$grad_dist(XX=XX)
          grad_samp <- MASS::mvrnorm(n=n, mu = grad_dist$mean[1,], Sigma = grad_dist$cov[1,,])
          grad_samp
          # gs2 <- apply(gs, 1, . %>% sum((.)^2))
          # c(mean(1/gs2), var(1/gs2))
        },
        grad_norm2_dist = function(XX) {
          # Calculate mean and var for squared norm of gradient
          # grad_dist <- gp$grad_dist(XX=XX) # Too slow because it does all
          d <- ncol(XX)
          nn <- nrow(XX)
          means <- numeric(nn)
          vars <- numeric(nn)
          for (i in 1:nn) {
            grad_dist_i <- self$grad_dist(XX=XX[i, , drop=FALSE])
            mean_i <- grad_dist_i$mean[1,]
            Sigma_i <- grad_dist_i$cov[1,,]
            SigmaInv_i <- solve(Sigma_i)
            # Using my own sqrt function since it is faster.
            # SigmaInvRoot_i <- expm::sqrtm(SigmaInv_i)
            SigmaInvRoot_i <- sqrt_matrix(mat=SigmaInv_i, symmetric = TRUE)
            eigen_i <- eigen(Sigma_i)
            P_i <- t(eigen_i$vectors)
            lambda_i <- eigen_i$values
            # testthat::expect_equal(t(P) %*% diag(eth$values) %*% (P), Sigma) # Should be equal
            b_i <- P_i %*% SigmaInvRoot_i %*% mean_i
            g2mean_i <- sum(b_i^2 * lambda_i)+d
            g2var_i <- 4*sum(b_i^2 * lambda_i^2)+2*d
            means[i] <- g2mean_i
            vars[i] <- g2var_i
          }
          data.frame(mean=means, var=vars)
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
