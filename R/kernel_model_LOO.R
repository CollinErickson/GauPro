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
    update = function (Xnew=NULL, Znew=NULL, Xall=NULL, Zall=NULL,
                       restarts = 5,
                       param_update = self$param.est, nug.update = self$nug.est, no_update=FALSE) {
      self$update_data(Xnew=Xnew, Znew=Znew, Xall=Xall, Zall=Zall) # Doesn't update Kinv, etc

      if (!no_update || (!param_update && !nug.update)) { # This option lets it skip parameter optimization entirely
        self$update_params(restarts=restarts, param_update=param_update,nug.update=nug.update)
      }

      self$update_K_and_estimates()


      # Do LOO stuff
      # browser()
      # Should I put this in an if use_LOO?
      # I don't want this not fit, then have user set use_LOO=T and have an error
      #   when it tries to predict with LOO.
      if (is.null(self$tmod)) {
        Zp <- self$pred_LOO(se.fit=TRUE)
        abs_t <- abs(Zp$t)
        self$tmod <- GauPro_kernel_model$new(X=self$X, Z=abs_t, kernel=GauPro::Matern52)
      } else {
        Zp <- self$pred_LOO(se.fit=TRUE)
        abs_t <- abs(Zp$t)
        self$tmod$update(Xall=self$X, Zall=abs_t, no_update=no_update)
      }

      invisible(self)
    },

    pred_one_matrix = function(XX, se.fit=F, covmat=F) {
      # input should already be check for matrix
      kxx <- self$kernel$k(XX) + self$nug
      kx.xx <- self$kernel$k(self$X, XX)
      mu_hat_matXX <- self$trend$Z(XX)

      mn <- pred_meanC_mumat(XX, kx.xx, self$mu_hatX, mu_hat_matXX, self$Kinv, self$Z)

      if (self$normalize) {mn <- mn * self$normalize_sd + self$normalize_mean}

      if (!se.fit & !covmat) {
        return(mn)
      }
      if (covmat) {
        stop("covmat not implemented for GauPro_kernel_model_LOO #68239")
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

      # browser()
      # Do LOO stuff here
      if (self$use_LOO) {
        loo.p <- self$tmod$predict(XX)
        se <- se * loo.p
        se <- pmax(se, 1e-8)
        s2 <- se ^ 2 #s2 * loo.p ^ 2
      }

      # se.fit but not covmat
      data.frame(mean=mn, s2=s2, se=se)
    }
  )
)
