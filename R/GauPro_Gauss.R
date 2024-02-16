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
#' @field corr Name of correlation
#' @field theta Correlation parameters
#' @field theta_length Length of theta
#' @field theta_map Map for theta
#' @field theta_short Short vector for theta
#' @field separable Are the dimensions separable?
#' @examples
#' n <- 12
#' x <- matrix(seq(0,1,length.out = n), ncol=1)
#' y <- sin(2*pi*x) + rnorm(n,0,1e-1)
#' gp <- GauPro_Gauss$new(X=x, Z=y, parallel=FALSE)
GauPro_Gauss <- R6::R6Class(
  classname = "GauPro_Gauss",
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
    #' @description Create GauPro object
    #' @param X Matrix whose rows are the input points
    #' @param Z Output points corresponding to X
    #' @param verbose Amount of stuff to print. 0 is little, 2 is a lot.
    #' @param separable Are dimensions separable?
    #' @param useC Should C code be used when possible? Should be faster.
    #' @param useGrad Should the gradient be used?
    #' @param parallel Should code be run in parallel? Make optimization
    #' faster but uses more computer resources.
    #' @param nug Value for the nugget. The starting value if estimating it.
    #' @param nug.min Minimum allowable value for the nugget.
    #' @param nug.est Should the nugget be estimated?
    #' @param param.est Should the kernel parameters be estimated?
    #' @param theta Correlation parameters
    #' @param theta_short Correlation parameters, not recommended
    #' @param theta_map Correlation parameters, not recommended
    #' @param ... Not used
    initialize = function(X, Z, verbose=0, separable=T, useC=F,useGrad=T,
                          parallel=FALSE,
                          nug=1e-6, nug.min=1e-8, nug.est=T,
                          param.est=T,
                          theta = NULL, theta_short = NULL, theta_map = NULL,
                          ...) {
      # This doesn't work inside R6 init as of now.
      # Deprecating, haven't used in years. Use kernel model instead.
      # Can't use deprecate_soft:
      #   https://github.com/r-lib/lifecycle/issues/149#issuecomment-1790953600
      lifecycle::deprecate_warn(
        when = "0.2.12",
        what = "GauPro_Gauss$new()",
        details = paste0("Please use GauPro::GauPro_kernel_model instead")
      )

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
    #' @description Correlation function
    #' @param x First point
    #' @param x2 Second point
    #' @param theta Correlation parameter
    corr_func = function(x, x2=NULL, theta=self$theta) {
      if (is.null(x2)) corr_gauss_matrix_symC(x, theta)
      else corr_gauss_matrixC(x, x2, theta)
    },
    #' @description Calculate deviance
    #' @param theta Correlation parameter
    deviance_theta = function (theta) {
      self$deviance(theta=theta, nug=self$nug)
    },
    #' @description Calculate deviance
    #' @param beta Correlation parameter on log scale
    deviance_theta_log = function (beta) {
      theta <- 10^beta
      self$deviance_theta(theta=theta)
    },
    #' @description Calculate deviance
    #' @param theta Correlation parameter
    #' @param nug Nugget
    deviance = function(theta=self$theta, nug=self$nug) {
      if (length(theta) < self$D) {theta = theta[self$theta_map]} # if not fully separable, map out to full theta
      #Gaussian_devianceC(theta, nug, self$X, self$Z) # adding a try in case C can't chol()
      try.dev.C <- try(Gaussian_devianceC(theta, nug, self$X, self$Z))
      if (inherits(try.dev.C, "try-error")) {return(Inf)}
      try.dev.C
    },
    #deviance_out = function(...){Gaussian_devianceC(...)},
    #' @description Calculate deviance gradient
    #' @param theta Correlation parameter
    #' @param nug Nugget
    #' @param joint Calculate over theta and nug at same time?
    #' @param overwhat Calculate over theta and nug at same time?
    deviance_grad = function(theta=NULL, nug=self$nug, joint=NULL,
                             overwhat=if (self$nug.est) "joint" else "theta") {
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
    #' @description Calculate deviance and gradient at same time
    #' @param theta Correlation parameter
    #' @param nug Nugget
    #' @param joint Calculate over theta and nug at same time?
    #' @param overwhat Calculate over theta and nug at same time?
    deviance_fngr = function (theta=NULL, nug=NULL,
                              overwhat=if (self$nug.est) "joint" else "theta") {
      if (length(theta) < self$D) {theta = theta[self$theta_map]} # if not fully separable, map out to full theta
      fngr <- Gaussian_deviance_fngrC(theta=theta, nug=nug, X=self$X, Z=self$Z, overwhat=overwhat)
      if (!self$separable & overwhat!="nug") { # Map down if theta not full dimensional
        gr <- fngr$gr
        tgr <- sapply(1:self$theta_length, function(ii){sum(gr[which(self$theta_map == ii)])})
        if (overwhat == "joint") { # If joint, add nugget to end and return
          return(list(fn=fngr$fn, gr=c(tgr, tail(gr,1))))
        }
        return(list(fn=fngr$fn, gr=tgr)) # else just return theta grad
      }
      fngr
    },
    #' @description Calculate deviance gradient
    #' @param beta Correlation parameter on log scale
    #' @param nug Nugget
    #' @param joint Calculate over theta and nug at same time?
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
    #' @description Calculate deviance on log scale
    #' @param beta Correlation parameter on log scale
    #' @param lognug Log of nugget
    #' @param joint Calculate over theta and nug at same time?
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
    #' @description Calculate deviance gradient on log scale
    #' @param beta Correlation parameter
    #' @param nug Nugget
    #' @param joint Calculate over theta and nug at same time?
    #' @param overwhat Calculate over theta and nug at same time?
    deviance_log_grad = function (beta=NULL, nug=self$nug, joint=NULL,
                                  overwhat=if (self$nug.est) "joint" else "theta") {
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
    #' @description Calculate deviance gradient on log scale
    #' @param beta Correlation parameter
    #' @param lognug Log of nugget
    #' @param joint Calculate over theta and nug at same time?
    #' @param overwhat Calculate over theta and nug at same time?
    deviance_log2_grad = function (beta=NULL, lognug=NULL, joint=NULL,
                                   overwhat=if (self$nug.est) "joint" else "theta") {
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
    #' @description Calculate deviance and gradient on log scale
    #' @param beta Correlation parameter
    #' @param lognug Log of nugget
    #' @param joint Calculate over theta and nug at same time?
    #' @param overwhat Calculate over theta and nug at same time?
    deviance_log2_fngr = function (beta=NULL, lognug=NULL, joint=NULL,
                                   overwhat=if (self$nug.est) "joint" else "theta") {
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
    #' @description Get optimization functions
    #' @param param_update Should the parameters be updated?
    #' @param nug.update Should the nugget be updated?
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
    #' @description Lower bound of params
    param_optim_lower = function() {# - lower bound of params
      rep(-5, self$theta_length)
    },
    #' @description Upper bound of params
    param_optim_upper = function() {# - upper bound of params
      rep(5, self$theta_length)
    },
    #' @description Start value of params for optim
    param_optim_start = function() {# - current param values on optim scale
      pmin(pmax(log(self$theta_short, 10), -5), 5)
    },
    #' @description Start value of params for optim
    param_optim_start0 = function () {# - some central param values that can be used for optimization restarts
      rep(0, self$theta_length)
    },
    #' @description Jitter value of params for optim
    #' @param param_value param value to add jitter to
    param_optim_jitter = function (param_value) { # only returns the jitter, not the param+jitter
      rnorm(self$theta_length,0,2)
    },
    #' @description Update value of params after optim
    #' @param restarts Number of restarts
    #' @param param_update Are the params being updated?
    #' @param nug.update Is the nugget being updated?
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
    #' @description Calculate the gradient
    #' @param XX Points to calculate grad at
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
    #' @description Calculate the gradient distribution
    #' @param XX Points to calculate grad at
    grad_dist = function (XX) {
      if (!is.matrix(XX)) {
        if (self$D == 1) XX <- matrix(XX, ncol=1)
        else if (length(XX) == self$D) XX <- matrix(XX, nrow=1)
        else stop('Predict input should be matrix')
      } else {
        if (ncol(XX) != self$D) {stop("Wrong dimension input")}
      }
      kx.xx <- self$corr_func(self$X, XX, theta=self$theta)

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
    #' @description Calculate the hessian
    #' @param XX Points to calculate grad at
    #' @param useC Should C code be used to speed up?
    hessian = function(XX, useC=self$useC) {
      if (!is.matrix(XX)) {
        if (self$D == 1) XX <- matrix(XX, ncol=1)
        else if (length(XX) == self$D) XX <- matrix(XX, nrow=1)
        else stop('Predict input should be matrix')
      } else {
        if (ncol(XX) != self$D) {stop("Wrong dimension input")}
      }
      hessian_func <- if (useC) {Gaussian_hessianC} else {Gaussian_hessianR}
      if (nrow(XX) == 1) {
        hess1 <- hessian_func(XX=XX, X=self$X, Z=self$Z, Kinv=self$Kinv, self$mu_hat, self$theta)
      } else { # more than one row
        hess1 <- lapply(1:nrow(XX), function(ii) {hessian_func(XX=XX[ii, ], X=self$X, Z=self$Z, Kinv=self$Kinv, self$mu_hat, self$theta)})
      }
      hess1
    },
    #' @description Print this object
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
