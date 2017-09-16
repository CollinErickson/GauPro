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
#' gp <- GauPro_Gauss_LOO$new(X=x, Z=y, parallel=FALSE)
GauPro_Gauss_LOO <- R6::R6Class(classname = "GauPro_Gauss_LOO",
  inherit = GauPro_Gauss,
  public = list(
    tmod = NULL,
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
      if (is.null(self$tmod)) {
        Zp <- self$pred_LOO(se.fit=TRUE)
        abs_t <- abs(Zp$t)
        self$tmod <- GauPro_Gauss$new(X=self$X, Z=abs_t)
      } else {
        Zp <- self$pred_LOO(se.fit=TRUE)
        abs_t <- abs(Zp$t)
        self$tmod$update(Xall=self$X, Zall=abs_t)
      }

      invisible(self)
    },

    pred_one_matrix = function(XX, se.fit=F, covmat=F) {
      # input should already be check for matrix
      kxx <- self$corr_func(XX) + self$nug
      kx.xx <- self$corr_func(self$X, XX)
      mn <- pred_meanC(XX, kx.xx, self$mu_hat, self$Kinv, self$Z)

      if (!se.fit & !covmat) {
        return(mn)
      }
      if (covmat) {
        stop("covmat not implemented for GauPro_Gauss_LOO #725934")
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

      # browser()
      # Do LOO stuff here
      loo.p <- self$tmod$predict(XX)
      se <- se * loo.p
      se <- pmax(se, 1e-8)
      s2 <- se ^ 2 #s2 * loo.p ^ 2


      # se.fit but not covmat
      data.frame(mean=mn, s2=s2, se=se)
    }
  )
)
