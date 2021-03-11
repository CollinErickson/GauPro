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
#' @field tmod A second GP model for the t-values of leave-one-out predictions
#' @field use_LOO Should the leave-one-out error corrections be used?
#' @examples
#' n <- 12
#' x <- matrix(seq(0,1,length.out = n), ncol=1)
#' y <- sin(2*pi*x) + rnorm(n,0,1e-1)
#' gp <- GauPro_kernel_model_LOO$new(X=x, Z=y, kernel=Gaussian)
#' y <- x^2 * sin(2*pi*x) + rnorm(n,0,1e-3)
#' gp <- GauPro_kernel_model_LOO$new(X=x, Z=y, kernel=Matern52)
#' y <- exp(-1.4*x)*cos(7*pi*x/2)
#' gp <- GauPro_kernel_model_LOO$new(X=x, Z=y, kernel=Matern52)
GauPro_kernel_model_LOO <- R6::R6Class(
  classname = "GauPro_kernel_model_LOO",
  inherit = GauPro_kernel_model,
  public = list(
    tmod = NULL, # A second GP model for the t-values of leave-one-out predictions
    use_LOO = TRUE, # Should predicted errors use leave-one-out correction?
    # Initialize: kernel, nug for LOO model, rest passed on to super$initialize
    #' @description Create a kernel model that uses a leave-one-out GP model
    #' to fix the standard error predictions.
    #' @param ... Passed to super$initialize.
    #' @param LOO_kernel The kernel that should be used for the
    #' leave-one-out model. Shouldn't be too smooth.
    #' @param LOO_options Options passed to the leave-one-out model.
    initialize = function (..., LOO_kernel, LOO_options=list()) {
      # Set tmod to be a function that will use input if given
      if (!missing(LOO_kernel) || !(length(LOO_options)==0)) {
        if (!missing(LOO_kernel)) {LOO_options$kernel <- LOO_kernel}
        else {LOO_options$kernel <- GauPro::Matern52}
        self$tmod <- function(X,Z) {
          LOO_options$X <- X
          LOO_options$Z <- Z
          do.call(GauPro_kernel_model$new, LOO_options)
        }
      }

      # Initialize rest with super
      super$initialize(...)
    },
    #' @description Update the model. Should only give in
    #' (Xnew and Znew) or (Xall and Zall).
    #' @param Xnew New X values to add.
    #' @param Znew New Z values to add.
    #' @param Xall All X values to be used. Will replace existing X.
    #' @param Zall All Z values to be used. Will replace existing Z.
    #' @param nug.update Is the nugget being updated?
    #' @param restarts Number of optimization restarts.
    #' @param param_update Are the parameters being updated?
    #' @param no_update Are no parameters being updated?
    update = function (Xnew=NULL, Znew=NULL, Xall=NULL, Zall=NULL,
                       restarts = 5,
                       param_update = self$param.est,
                       nug.update = self$nug.est, no_update=FALSE) {
      self$update_data(Xnew=Xnew, Znew=Znew, Xall=Xall, Zall=Zall) # Doesn't update Kinv, etc

      if (!no_update || (!param_update && !nug.update)) { # This option lets it skip parameter optimization entirely
        self$update_params(restarts=restarts, param_update=param_update,nug.update=nug.update)
      }

      self$update_K_and_estimates()


      # Do LOO stuff
      # Should I put this in an if use_LOO?
      #   I don't want this not fit, then have user set use_LOO=T and have an error
      #   when it tries to predict with LOO.
      if (is.null(self$tmod)) {
        Zp <- self$pred_LOO(se.fit=TRUE)
        abs_t <- abs(Zp$t)
        self$tmod <- GauPro_kernel_model$new(X=self$X, Z=abs_t, kernel=GauPro::Matern52)
      } else if (is.function(self$tmod)) { # Premade function, only call with X and Z
        Zp <- self$pred_LOO(se.fit=TRUE)
        abs_t <- abs(Zp$t)
        self$tmod <- self$tmod(X=self$X, Z=abs_t)
      } else {
        Zp <- self$pred_LOO(se.fit=TRUE)
        abs_t <- abs(Zp$t)
        self$tmod$update(Xall=self$X, Zall=abs_t, no_update=no_update)
      }

      invisible(self)
    },
    #' @description Predict for a matrix of points
    #' @param XX points to predict at
    #' @param se.fit Should standard error be returned?
    #' @param covmat Should covariance matrix be returned?
    #' @param return_df When returning se.fit, should it be returned in
    #' a data frame?
    pred_one_matrix = function(XX, se.fit=F, covmat=F, return_df=FALSE) {
      # input should already be check for matrix
      # kxx <- self$kernel$k(XX) + diag(self$nug * self$s2_hat, nrow(XX))
      kx.xx <- self$kernel$k(self$X, XX)
      mu_hat_matXX <- self$trend$Z(XX)

      mn <- pred_meanC_mumat(XX, kx.xx, self$mu_hatX, mu_hat_matXX, self$Kinv, self$Z)

      if (self$normalize) {mn <- mn * self$normalize_sd + self$normalize_mean}

      if (!se.fit & !covmat) {
        return(mn)
      }
      if (covmat && self$use_LOO) {
        stop("covmat not implemented for GauPro_kernel_model_LOO #68239")
        # new for kernel
        kxx <- self$kernel$k(XX) + diag(self$nug * self$s2_hat, nrow(XX))
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
      # covmatdat <- kxx - t(kx.xx) %*% self$Kinv %*% kx.xx
      # s2 <- diag(covmatdat)
      # Speed up, see kernel_model.R for more details
      diag.kxx <- self$nug * self$s2_hat + rep(self$s2_hat, nrow(XX))
      s2 <- diag.kxx - colSums( (kx.xx) * (self$Kinv %*% kx.xx))

      if (self$normalize) {
        s2 <- s2 * self$normalize_sd ^ 2
      }

      # s2 <- pred_var(XX, kxx, kx.xx, self$s2_hat, self$Kinv, self$Z)
      se <- rep(0, length(mn)) # NEG VARS will be 0 for se, NOT SURE I WANT THIS
      se[s2>=0] <- sqrt(s2[s2>=0])

      # browser()
      # Do LOO stuff here
      if (self$use_LOO) {
        loo.p <- self$tmod$predict(XX)
        se <- se * loo.p
        se <- pmax(se, 1e-8)
        s2 <- se ^ 2 #s2 * loo.p ^ 2
      }

      # se.fit but not covmat
      if (return_df) {
        data.frame(mean=mn, s2=s2, se=se) # data.frame is really slow compared to cbind or list
      } else {
        list(mean=mn, s2=s2, se=se)
      }
    }
  )
)
