#' GauPro model that uses kernels
#'
#' Class providing object with methods for fitting a GP model.
#' Allows for different kernel and trend functions to be used.
#'
#' @docType class
#' @importFrom R6 R6Class
#' @importFrom stats model.frame
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
#' gp$plot1D()
#' gp$cool1Dplot()
#'
#' n <- 200
#' d <- 7
#' x <- matrix(runif(n*d), ncol=d)
#' f <- function(x) {x[1]*x[2] + cos(x[3]) + x[4]^2}
#' y <- apply(x, 1, f)
#' gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian)
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
# @field corr Type of correlation function
#' @field nug.min Minimum value of nugget
#' @field nug.max Maximum value of the nugget.
#' @field nug.est Should the nugget be estimated?
#' @field nug Value of the nugget, is estimated unless told otherwise
# @field separable Are the dimensions separable?
#' @field param.est Should the kernel parameters be estimated?
#' @field verbose 0 means nothing printed, 1 prints some, 2 prints most.
#' @field useGrad Should grad be used?
#' @field useC Should C code be used?
#' @field parallel Should the code be run in parallel?
#' @field parallel_cores How many cores are there? By default it detects.
#' @field kernel The kernel to determine the correlations.
#' @field trend The trend.
#' @field mu_hatX Predicted trend value for each point in X.
#' @field s2_hat Variance parameter estimate
#' @field K Covariance matrix
#' @field Kchol Cholesky factorization of K
#' @field Kinv Inverse of K
#' @field Kinv_Z_minus_mu_hatX K inverse times Z minus the predicted
#' trend at X.
#' @field restarts Number of optimization restarts to do when updating.
#' @field normalize Should the inputs be normalized?
#' @field normalize_mean If using normalize, the mean of each column.
#' @field normalize_sd If using normalize, the standard
#' deviation of each column.
#' @field optimizer What algorithm should be used to optimize the
#' parameters.
#' @field track_optim Should it track the parameters evaluated
#' while optimizing?
#' @field track_optim_inputs If track_optim is TRUE,
#' this will keep a list of parameters evaluated.
#' View them with plot_track_optim.
#' @field track_optim_dev If track_optim is TRUE,
#' this will keep a vector of the deviance values calculated
#' while optimizing parameters.
#' View them with plot_track_optim.
#' @field formula Formula
#' @field convert_formula_data List for storing data to convert data
#' using the formula
#' @section Methods:
#' \describe{
#'   \item{\code{new(X, Z, corr="Gauss", verbose=0, separable=T, useC=F,
#'                   useGrad=T,
#'          parallel=T, nug.est=T, ...)}}{
#'          This method is used to create object of this
#'          class with \code{X} and \code{Z} as the data.}
#'
#'   \item{\code{update(Xnew=NULL, Znew=NULL, Xall=NULL, Zall=NULL,
#' restarts = 0,
#' param_update = T, nug.update = self$nug.est)}}{This method updates the
#' model, adding new data if given, then running optimization again.}
#'   }
# GauPro_kernel_model ----
GauPro_kernel_model <- R6::R6Class(
  classname = "GauPro",
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
    track_optim = NULL,
    track_optim_inputs = list(),
    track_optim_dev = numeric(0),
    formula = NULL,
    convert_formula_data = NULL,
    #deviance_out = NULL, #(theta, nug)
    #deviance_grad_out = NULL, #(theta, nug, overwhat)
    #deviance_fngr_out = NULL,
    #' @description Create kernel_model object
    #' @param X Matrix whose rows are the input points
    #' @param Z Output points corresponding to X
    #' @param kernel The kernel to use. E.g., Gaussian$new().
    #' @param trend Trend to use. E.g., trend_constant$new().
    #' @param verbose Amount of stuff to print. 0 is little, 2 is a lot.
    #' @param useC Should C code be used when possible? Should be faster.
    #' @param useGrad Should the gradient be used?
    #' @param parallel Should code be run in parallel? Make optimization
    #' faster but uses more computer resources.
    #' @param parallel_cores When using parallel, how many cores should
    #' be used?
    #' @param nug Value for the nugget. The starting value if estimating it.
    #' @param nug.min Minimum allowable value for the nugget.
    #' @param nug.max Maximum allowable value for the nugget.
    #' @param nug.est Should the nugget be estimated?
    #' @param param.est Should the kernel parameters be estimated?
    #' @param restarts How many optimization restarts should be used when
    #' estimating parameters?
    #' @param normalize Should the data be normalized?
    #' @param optimizer What algorithm should be used to optimize the
    #' parameters.
    #' @param track_optim Should it track the parameters evaluated
    #' while optimizing?
    #' @param ... Not used
    initialize = function(X, Z,
                          kernel, trend,
                          verbose=0, useC=TRUE, useGrad=TRUE,
                          parallel=FALSE, parallel_cores="detect",
                          nug=1e-6, nug.min=1e-8, nug.max=Inf, nug.est=TRUE,
                          param.est = TRUE, restarts = 0,
                          normalize = FALSE, optimizer="L-BFGS-B",
                          track_optim=FALSE,
                          ...) {
      #self$initialize_GauPr(X=X,Z=Z,verbose=verbose,useC=useC,
      #                      useGrad=useGrad,
      #                      parallel=parallel, nug.est=nug.est)
      # if (missing(X)) {
      #   stop("You must give X to GauPro_kernel_model")
      # }
      # if (missing(Z)) {
      #   stop("You must give Z to GauPro_kernel_model")
      # }
      if ((!missing(X) && is.formula(X)) || (!missing(Z) && is.formula(Z))) {
        if (!missing(X) && is.formula(X)) {
          formula <- X
          # stopifnot(is.data.frame(Z))
          if (!missing(Z) && is.data.frame(Z)) {
            data <- Z
          } else {
            data <- NULL
          }
        }
        if (!missing(Z) && is.formula(Z)) {
          formula <- Z
          # stopifnot(is.data.frame(X))
          if (!missing(X) && is.data.frame(X)) {
            data <- X
          } else {
            data <- NULL
          }
        }
        modfr <- model.frame(formula = formula, data = data)
        Z <- modfr[,1]
        Xdf <- modfr[,2:ncol(modfr), drop=FALSE]
        convert_formula_data <- list(factors=list(),
                                     chars=list())
        # Convert factor columns to integer
        for (i in 1:ncol(Xdf)) {
          if (is.factor(Xdf[, i])) {
            convert_formula_data$factors[[
              length(convert_formula_data$factors)+1
            ]] <- list(index=i,
                       levels=levels(Xdf[[i]]))
            Xdf[[i]] <- as.integer(Xdf[[i]])
          }
        }
        # Convert char columns to integer
        for (i in 1:ncol(Xdf)) {
          if (is.character(Xdf[, i])) {
            convert_formula_data$chars[[
              length(convert_formula_data$chars)+1
            ]] <- list(index=i,
                       vals=sort(unique(Xdf[[i]])))
            Xdf[[i]] <- sapply(Xdf[[i]],
                               function(x) {
                                 which(x==convert_formula_data$chars[[
                                   length(convert_formula_data$chars)
                                 ]]$vals)
                               })
          }
        }
        self$formula <- formula
        # self$data <- data
        self$convert_formula_data <- convert_formula_data
        X <- as.matrix(Xdf)
      }

      if (missing(X) || is.null(X)) {
        stop("You must give X to GauPro_kernel_model")
      }
      if (missing(Z) || is.null(Z)) {
        stop("You must give Z to GauPro_kernel_model")
      }

      if (is.data.frame(X)) {
        X <- as.matrix(X)
      }
      stopifnot(is.numeric(X))
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
      if ("R6ClassGenerator" %in% class(kernel)) {
        # Let generator be given so D can be set auto
        self$kernel <- kernel$new(D=self$D)
      } else if ("GauPro_kernel" %in% class(kernel)) {
        # Otherwise it should already be a kernel
        self$kernel <- kernel
      } else if(is.character(kernel) && length(kernel)==1) {
        kernel <- tolower(kernel)
        if (kernel %in% c("gaussian", "gauss")) {
          kernel <- Gaussian$new(D=self$D, useC=useC)
        } else if (kernel %in% c("matern32", "m32", "matern3/2",
                                 "matern3_2")) {
          kernel <- Matern32$new(D=self$D)
        } else if (kernel %in% c("matern52", "m52", "matern5/2",
                                 "matern5_2")) {
          kernel <- Matern52$new(D=self$D)
        } else if (kernel %in% c("exp", "exponential",
                                 "m12", "matern12",
                                 "matern1/2", "matern1_2")) {
          kernel <- Exponential$new(D=self$D)
        } else if (kernel %in% c("ratquad", "rationalquadratic", "rq")) {
          kernel <- RatQuad$new(D=self$D)
        } else if (kernel %in% c("powerexponential", "powexp", "pe")) {
          kernel <- PowerExp$new(D=self$D)
        } else {
          stop(paste0("Kernel given to GauPro_kernel_model (",
                      kernel, ") is not valid. ",
                      "Consider using Gaussian or Matern52."))
        }
        self$kernel <- kernel
      } else {
        stop(paste0("Kernel given to GauPro_kernel_model is not valid. ",
                    "Consider using Gaussian or Matern52."))
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
      stopifnot(length(track_optim) == 1, is.logical(track_optim))
      self$track_optim <- track_optim

      self$update_K_and_estimates() # Need to get mu_hat before starting
      # self$mu_hat <- mean(Z)
      self$fit()
      invisible(self)
    },
    # initialize_GauPr = function() {
    # },
    #' @description  Fit model
    #' @param X Inputs
    #' @param Z Outputs
    fit = function(X, Z) {
      self$update()
    },
    #' @description Update covariance matrix and estimates
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
        cat("Increasing nugget to get invertibility from ", oldnug, ' to ', self$nug, "\n")
      }
      self$Kinv <- chol2inv(self$Kchol)
      # self$mu_hat <- sum(self$Kinv %*% self$Z) / sum(self$Kinv)
      self$mu_hatX <- self$trend$Z(X=self$X)
      self$Kinv_Z_minus_mu_hatX <- c(self$Kinv %*% (self$Z - self$mu_hatX))
      # self$s2_hat <- c(t(self$Z - self$mu_hat) %*% self$Kinv %*%
      #                               (self$Z - self$mu_hat) / self$N)
      self$s2_hat <- self$kernel$s2
    },
    #' @description Predict for a matrix of points
    #' @param XX points to predict at
    #' @param se.fit Should standard error be returned?
    #' @param covmat Should covariance matrix be returned?
    #' @param split_speed Should the matrix be split for faster predictions?
    #' @param mean_dist Should the error be for the distribution of the mean?
    predict = function(XX, se.fit=F, covmat=F, split_speed=F, mean_dist=FALSE) {
      self$pred(XX=XX, se.fit=se.fit, covmat=covmat,
                split_speed=split_speed, mean_dist=mean_dist)
    },
    #' @description Predict for a matrix of points
    #' @param XX points to predict at
    #' @param se.fit Should standard error be returned?
    #' @param covmat Should covariance matrix be returned?
    #' @param split_speed Should the matrix be split for faster predictions?
    #' @param mean_dist Should the error be for the distribution of the mean?
    pred = function(XX, se.fit=F, covmat=F, split_speed=F, mean_dist=FALSE) {
      if (!is.null(self$formula)) {
        XX <- convert_X_with_formula(XX, self$convert_formula_data,
                                     self$formula)
      }
      if (is.data.frame(XX)) {
        XX <- as.matrix(XX)
      }
      if (is.matrix(XX)) {
        stopifnot(is.numeric(XX))
      } else {
        if (is.numeric(XX)) {
          if (self$D == 1) XX <- matrix(XX, ncol=1)
          else if (length(XX) == self$D) XX <- matrix(XX, nrow=1)
          else stop(paste0('Predict input should be matrix with ', self$D,
                           ' columns or vector of length ', self$D))
        } else {
          stop(paste("Bad type of XX given to pred"))
        }
      }

      stopifnot(ncol(XX) == self$D)
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
                                        covmat=covmat, mean_dist=mean_dist)
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
                                      covmat=covmat, return_df=TRUE,
                                      mean_dist=mean_dist)
        return(pred1)
      }
    },
    #' @description Predict for a matrix of points
    #' @param XX points to predict at
    #' @param se.fit Should standard error be returned?
    #' @param covmat Should covariance matrix be returned?
    #' @param return_df When returning se.fit, should it be returned in
    #' a data frame?
    #' @param mean_dist Should the error be for the distribution of the mean?
    pred_one_matrix = function(XX, se.fit=F, covmat=F, return_df=FALSE,
                               mean_dist=FALSE) {
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
      # It's supposed to return a vector, but it's a matrix
      mn <- mn[, 1]

      if (self$normalize) {
        mn <- mn * self$normalize_sd + self$normalize_mean
      }

      if (!se.fit & !covmat) {
        return(mn)
      }
      if (covmat) {
        # new for kernel
        # kxx <- self$kernel$k(XX) + diag(self$nug * self$s2_hat, nrow(XX))
        kxx <- self$kernel$k(XX)
        # The mean doesn't get the nugget added
        if (!mean_dist) {
          kxx <- kxx + diag(self$nug * self$s2_hat, nrow(XX))
        }
        covmatdat <- kxx - t(kx.xx) %*% self$Kinv %*% kx.xx

        if (self$normalize) {
          covmatdat <- covmatdat * self$normalize_sd ^ 2
        }

        # #covmatdat <- self$pred_var(XX, kxx=kxx, kx.xx=kx.xx, covmat=T)
        # covmatdat <- pred_cov(XX, kxx, kx.xx, self$s2_hat, self$Kinv,
        #                       self$Z)
        s2 <- diag(covmatdat)
        # se <- rep(1e-8, length(mn)) # NEG VARS will be 0 for se,
        # #  NOT SURE I WANT THIS
        # se[s2>=0] <- sqrt(s2[s2>=0])

        if (any(s2 < 0)) {
          min_s2 <- .Machine$double.eps
          warning(paste0("Negative s2 predictions are being set to ",
                         min_s2, " (", sum(s2<0)," values).",
                         " covmat is not being altered."))
          s2 <- pmax(s2, min_s2)
        }
        se <- sqrt(s2)
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
      # diag.kxx <- self$nug * self$s2_hat + rep(self$s2_hat, nrow(XX))
      diag.kxx <- rep(self$s2_hat, nrow(XX))
      if (!mean_dist) {
        diag.kxx <- diag.kxx + self$nug * self$s2_hat
      }
      s2 <- diag.kxx - colSums( (kx.xx) * (self$Kinv %*% kx.xx))

      if (self$normalize) {
        s2 <- s2 * self$normalize_sd ^ 2
      }

      # # s2 <- pred_var(XX, kxx, kx.xx, self$s2_hat, self$Kinv, self$Z)
      # se <- rep(0, length(mn)) # NEG VARS will be 0 for se,
      # #   NOT SURE I WANT THIS
      # se[s2>=0] <- sqrt(s2[s2>=0])
      if (any(s2 < 0)) {
        min_s2 <- .Machine$double.eps
        warning(paste0("Negative s2 predictions are being set to ",
                       min_s2, " (", sum(s2<0)," values)"))
        s2 <- pmax(s2, min_s2)
      }
      se <- sqrt(s2)

      # se.fit but not covmat
      if (return_df) {
        # data.frame is really slow compared to cbind or list
        data.frame(mean=mn, s2=s2, se=se)
      } else {
        list(mean=mn, s2=s2, se=se)
      }
    },
    #' @description Predict mean
    #' @param XX points to predict at
    #' @param kx.xx Covariance of X with XX
    pred_mean = function(XX, kx.xx) { # 2-8x faster to use pred_meanC
      # c(self$mu_hat + t(kx.xx) %*% self$Kinv %*% (self$Z - self$mu_hat))
      # mu_hat_matX <- self$trend$Z(self$X)
      mu_hat_matXX <- self$trend$Z(XX)
      c(mu_hat_matXX + t(kx.xx) %*% self$Kinv %*% (self$Z - self$mu_hatX))
    },
    #' @description Predict mean using C
    #' @param XX points to predict at
    #' @param kx.xx Covariance of X with XX
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
    #' @description Predict variance
    #' @param XX points to predict at
    #' @param kxx Covariance of XX with itself
    #' @param kx.xx Covariance of X with XX
    #' @param covmat Should the covariance matrix be returned?
    pred_var = function(XX, kxx, kx.xx, covmat=F) {
      # 2-4x faster to use C functions pred_var and pred_cov
      self$s2_hat * diag(kxx - t(kx.xx) %*% self$Kinv %*% kx.xx)
    },
    #' @description leave one out predictions
    #' @param se.fit Should standard errors be included?
    pred_LOO = function(se.fit=FALSE) {
      # Predict LOO (leave-one-out) on data used to fit model
      # See vignette for explanation of equations
      # If se.fit==T, then calculate the LOO se and the
      #    corresponding t score
      Z_LOO <- numeric(self$N)
      if (se.fit) {Z_LOO_s2 <- numeric(self$N)}
      Z_trend <- self$trend$Z(self$X)
      for (i in 1:self$N) {
        E <- self$Kinv[-i, -i] # Kinv without i
        b <- self$K[    i, -i] # K    between i and rest
        g <- self$Kinv[ i, -i] # Kinv between i and rest
        # Kinv for K if i wasn't in K
        Ainv <- E + E %*% b %*% g / (1-sum(g*b))
        Zi_LOO <- Z_trend[i] + c(b %*% Ainv %*% (self$Z[-i] - Z_trend[-i]))
        Z_LOO[i] <- Zi_LOO
        if (se.fit) {
          Zi_LOO_s2 <- self$K[i,i] - c(b %*% Ainv %*% b)
          # Have trouble when s2 < 0, set to small number
          Zi_LOO_s2 <- max(Zi_LOO_s2, 1e-16)
          Z_LOO_s2[i] <- Zi_LOO_s2
        }
      }
      if (se.fit) { # Return df with se and t if se.fit
        Z_LOO_se <- sqrt(Z_LOO_s2)
        t_LOO <- (self$Z - Z_LOO) / Z_LOO_se
        data.frame(fit=Z_LOO, se.fit=Z_LOO_se, t=t_LOO)
      } else { # Else just mean LOO
        Z_LOO
      }
    },
    #' @description Predict variance after adding points
    #' @param add_points Points to add
    #' @param pred_points Points to predict at
    pred_var_after_adding_points = function(add_points, pred_points) {
      # Calculate pred_var at pred_points after add_points
      #  have been added to the design self$X
      # S is add points
      # G <- solve(self$pred(add_points, covmat = TRUE)$cov)
      # FF <- -self$Kinv %*% 1
      if (!is.matrix(add_points)) {
        if (length(add_points) != self$D) {
          stop("add_points must be matrix or of length D")
        }
        else {add_points <- matrix(add_points, nrow=1)}
      } else if (ncol(add_points) != self$D) {
        stop("add_points must have dimension D")
      }
      C_S <- self$kernel$k(add_points)
      C_S <- C_S + self$s2_hat * diag(self$nug, nrow(C_S)) # Add nugget
      C_XS <- self$kernel$k(self$X, add_points)
      C_X_inv_C_XS <- self$Kinv %*% C_XS
      G <- solve(C_S - t(C_XS) %*% C_X_inv_C_XS)
      FF <- - C_X_inv_C_XS %*% G
      E <- self$Kinv - C_X_inv_C_XS %*% t(FF)

      # Speed this up a lot by avoiding apply and doing all at once
      # Assume single point cov is s2(1+nug)
      C_a <- self$s2_hat * (1 + self$nug)
      C_Xa <- self$kernel$k(self$X, pred_points) # length n vec, not matrix
      C_Sa <- self$kernel$k(add_points, pred_points)
      C_a - (colSums(C_Xa * (E %*% C_Xa)) +
               2 * colSums(C_Xa * (FF %*% C_Sa)) +
               colSums(C_Sa * (G %*% C_Sa)))
    },
    #' @description Predict variance reductions after adding each point separately
    #' @param add_points Points to add
    #' @param pred_points Points to predict at
    pred_var_after_adding_points_sep = function(add_points, pred_points) {
      # Calculate pred_var at pred_points after each add_points
      #  has individually (separately) been added to the design self$X
      # A vectorized version of pred_var_after_adding_points_sep
      # S is add points, a is pred_points in variables below
      # Output is matrix of size nrow(pred_points) by nrow(add_points)
      #  where (i,j) element is predictive variance at pred_points[i,]
      #  after add_points[j,] has been added to current design
      # Equations below, esp with sweep and colSums are confusing
      #  but work out. Make it fast and vectorized.
      # Check against pred_var_after_adding_points.
      if (!is.matrix(add_points)) {
        if (length(add_points) != self$D) {
          stop("add_points must be matrix or of length D")
        }
        else {add_points <- matrix(add_points, nrow=1)}
      } else if (ncol(add_points) != self$D) {
        stop("add_points must have dimension D")
      }
      C_S <- self$s2_hat * (1+self$nug)
      C_XS <- self$kernel$k(self$X, add_points)
      C_X_inv_C_XS <- self$Kinv %*% C_XS
      G <- 1 / (C_S - colSums(C_XS * C_X_inv_C_XS))

      # Speed this up a lot by avoiding apply and doing all at once
      # Assume single point cov is s2(1+nug)
      C_a <- self$s2_hat * (1 + self$nug)
      C_Xa <- self$kernel$k(self$X, pred_points) # matrix
      C_Sa <- self$kernel$k(add_points, pred_points) # matrix
      t1a <- colSums(C_Xa * (self$Kinv %*% C_Xa))
      t1 <- sweep(sweep((t(C_Xa) %*% C_X_inv_C_XS)^2, 2, G, `*`),
                  1, t1a, `+`)
      t2 <- -2*sweep((t(C_X_inv_C_XS) %*% C_Xa) * C_Sa, 1, G, `*`)
      t3 <- sweep((C_Sa)^2, 1, G, `*`)
      return(C_a - (t1 + t(t2 + t3)))
    },
    #' @description Predict variance reduction for a single point
    #' @param add_point Point to add
    #' @param pred_points Points to predict at
    pred_var_reduction = function(add_point, pred_points) {
      # Calculate pred_var at pred_points after add_point
      #  have been added to the design self$X
      # S is add point
      if (!is.vector(add_point) || length(add_point)!=self$D) {
        stop("add_point must be vector of length D")
      }
      C_S <- self$s2_hat * (1 + self$nug) # Assumes correlation structure
      C_XS <- self$kernel$k(self$X, add_point)
      C_X_inv_C_XS <- as.vector(self$Kinv %*% C_XS)
      G <- 1 / c(C_S - t(C_XS) %*% C_X_inv_C_XS)

      # Assume single point cov is s2(1+nug)
      # C_a <- self$s2_hat * (1 + self$nug)
      # pred_var_a_func <- function(a) {
      #   C_Xa <- self$kernel$k(self$X, a) # length n vector, not matrix
      #   C_Sa <- self$kernel$k(add_point, a)
      #   (sum(C_Xa * C_X_inv_C_XS) - C_Sa) ^ 2 * G
      # }
      # if (is.matrix(pred_points)) {
      #   if (method1) { # Slow way
      #     prds <- apply(pred_points, 1, pred_var_a_func)
      #   } else {
      # Speeding it up by getting all at once instead of by row
      C_aX <- self$kernel$k(pred_points, self$X) # len n vec, not mat
      C_aS <- self$kernel$k(pred_points, add_point)
      # (sum(C_Xa * C_X_inv_C_XS) - C_Sa) ^ 2 * G
      prds <- (c(C_aX %*% C_X_inv_C_XS) - C_aS) ^ 2 * G
      # }
      # }
      # else {prds <- pred_var_a_func(pred_points)}
      prds
    },
    #' @description Predict variance reductions
    #' @param add_points Points to add
    #' @param pred_points Points to predict at
    pred_var_reductions = function(add_points, pred_points) {
      # Calculate pred_var at pred_points after each of add_points
      #  has been added to the design self$X separately.
      # This is a vectorized version of pred_var_reduction,
      #  to consider all of add_points added together
      #  use pred_var_after_adding_points
      # S is add points, a is pred points
      if (!is.matrix(add_points) || ncol(add_points) != self$D) {
        stop("add_points must be a matrix with D columns")
      }
      # C_S <- self$s2_hat * diag(1 + self$nug, nrow(add_points))
      # Assumes correlation structure
      C_XS <- self$kernel$k(self$X, add_points)
      C_X_inv_C_XS <- self$Kinv %*% C_XS
      # G <- 1 / c(C_S - t(C_XS) %*% C_X_inv_C_XS)
      G <- 1 / (self$s2_hat*(1+self$nug) - colSums(C_XS * C_X_inv_C_XS))

      C_aX <- self$kernel$k(pred_points, self$X) # now a matrix
      C_aS <- self$kernel$k(pred_points, add_points)
      # (sum(C_Xa * C_X_inv_C_XS) - C_Sa) ^ 2 * G
      prds <- sweep(((C_aX %*% C_X_inv_C_XS) - C_aS) ^ 2, 2, G, `*`)
      prds
    },
    #' @description Plot the object
    #' @param ... Parameters passed to cool1Dplot(), plot2D(), or plotmarginal()
    plot = function(...) {
      if (self$D == 1) {
        self$cool1Dplot(...)
      } else if (self$D == 2) {
        self$plot2D(...)
      } else {
        # stop("No plot method for higher than 2 dimension")
        self$plotmarginal(...)
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
      px <- self$pred(x, covmat = T, mean_dist=TRUE)
      # px$cov <- self$kernel$k(matrix(x,ncol=1))
      # n2 <- 20
      Sigma.try <- try(newy <- MASS::mvrnorm(n=n2, mu=px$mean,
                                             Sigma=px$cov),
                       silent = TRUE)
      nug_Sig <- 1e-8 # self$nug, now use 1e-8 since self$nug is excluded in pred.
      while (inherits(Sigma.try, "try-error")) {
        message(paste0("Adding nugget to cool1Dplot: ", nug_Sig))
        Sigma.try <- try(
          newy <- MASS::mvrnorm(n=n2, mu=px$mean,
                                Sigma=px$cov + diag(nug_Sig, nrow(px$cov))),
          silent = TRUE)
        # if (inherits(Sigma.try2, "try-error")) {
        #   stop("Can't do cool1Dplot")
        # }
        nug_Sig <- 2*nug_Sig
      }
      if (n2==1) { # Avoid error when n2=1
        newy <- matrix(newy, nrow=1)
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
             if (self$normalize) {
               self$Z * self$normalize_sd + self$normalize_mean
             } else {self$Z},
             pch=19, col=1, cex=2)
    },
    #' @description Make 1D plot
    #' @param n2 Number of things to plot
    #' @param nn Number of things to plot
    #' @param col2 Color of the prediction interval
    #' @param col3 Color of the interval for the mean
    #' @param ylab y label
    #' @param xlab x label
    #' @param xmin xmin
    #' @param xmax xmax
    #' @param ymax ymax
    #' @param ymin ymin
    plot1D = function(n2=20, nn=201, col2=2, col3=3, #"gray",
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
      pxmean <- self$pred(x, se=T, mean_dist=T)
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
           xlab=xlab, ylab=ylab,
           main=paste0("Predicted output (95% interval for mean is green,",
                       " 95% interval for sample is red)"))
      points(x, px$mean-2*px$se, type='l', col=col2, lwd=2)
      points(x, pxmean$mean+2*pxmean$se, type='l', col=col3, lwd=2)
      points(x, pxmean$mean-2*pxmean$se, type='l', col=col3, lwd=2)
      points(x,px$me, type='l', lwd=4)
      points(self$X,
             if (self$normalize) {
               self$Z * self$normalize_sd + self$normalize_mean
             } else {self$Z},
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
      ContourFunctions::cf_func(self$predict, batchmax=Inf,
                                xlim=c(xmin, xmax),
                                ylim=c(ymin, ymax),
                                pts=self$X,
                                gg=TRUE)
    },
    #' @description Plot marginal. For each input, hold all others at a constant
    #' value and adjust it along it's range to see how the prediction changes.
    #' @param pt What point to use as a center. All values except the one being
    #' plotted are held constant at this value.
    plotmarginal = function(pt=colMeans(self$X)) {
      # pt <- colMeans(self$X)
      # pt
      factorinfo <- find_kernel_factor_dims(self$kernel)
      if (length(factorinfo > 0)) {
        factorindexes <- factorinfo[2*(1:(length(factorinfo)/2))-1]
        factornlevels <- factorinfo[2*(1:(length(factorinfo)/2))]
        for (i in 1:length(factorindexes)) {
          if (!(pt[factorindexes[i]] %in% 1:factornlevels[i])) {
            pt[factorindexes[i]] <- sample(1:factornlevels[i], 1)
          }
        }
      } else {
        factorindexes <- c()
      }
      pts <- NULL
      for (i in 1:ncol(self$X)) {
        if (i %in% factorindexes) {
          ind_i <- which(factorindexes == i)
          xseq <- 1:(factorinfo[2*ind_i])
        } else {
          xseq <- seq(min(self$X[,i]), max(self$X[,i]), l=51)
        }
        Xmat <- matrix(pt, byrow=T, ncol=length(pt), nrow=length(xseq))
        Xmat[, i] <- xseq
        pX <- self$pred(Xmat, se.fit = T)
        pts <- rbind(pts,
                     cbind(pred=pX$mean, predse=pX$se, xi=xseq, i=i))
      }
      pts2 <- as.data.frame(pts)
      # pts2 %>%
      #   mutate(predupper=pred+2*predse,
      #          predlower=pred-2*predse)
      pts2$predupper <- pts2$pred + 2*pts2$predse
      pts2$predlower <- pts2$pred - 2*pts2$predse
      ggplot2::ggplot(data=pts2, ggplot2::aes(xi, pred)) +
        ggplot2::facet_wrap(.~i, scales = "free_x") +
        ggplot2::geom_line(ggplot2::aes(y=predupper), color="green") +
        ggplot2::geom_line(ggplot2::aes(y=predlower), color="green") +
        ggplot2::geom_line(size=1) +
        ggplot2::ylab("Predicted Z (95% interval)") +
        ggplot2::xlab("x along dimension i")
    },
    #' @description Plot marginal prediction for random sample of inputs
    #' @param n Number of random points to evaluate
    plotmarginalrandom = function(n=151) {
      # Plot marginal random averages
      # Get random matrix, scale to proper lower/upper
      X <- lhs::randomLHS(n=n, k=ncol(self$X))
      X2 <- sweep(X, 2, apply(self$X, 2, max) - apply(self$X, 2, min), "*")
      X3 <- sweep(X2, 2, apply(self$X, 2, min), "+")
      factorinfo <- find_kernel_factor_dims(self$kernel)
      if (length(factorinfo) > 0) {
        for (i in 1:(length(factorinfo)/2)) {
          X3[, factorinfo[i*2-1]] <- sample(1:factorinfo[i*2], size=n, replace=T)
        }
      }
      X3pred <- self$pred(X3, se.fit = T)
      X3pred$irow <- 1:nrow(X3pred)
      # head(X3pred)
      X4 <- dplyr::inner_join(
        X3pred,
        tidyr::pivot_longer(cbind(as.data.frame(X3),irow=1:nrow(X)), cols=1:ncol(self$X)),
        "irow")
      # head(X4)
      X4$upper <- X4$mean + 2*X4$se
      X4$lower <- X4$mean - 2*X4$se
      ggplot2::ggplot(X4, ggplot2::aes(value, mean)) +
        ggplot2::facet_wrap(.~name, scales="free_x") +
        # geom_point(aes(y=upper), color="green") +
        ggplot2::geom_segment(ggplot2::aes(y=upper, yend=lower, xend=value),
                              color="green", size=2) +
        ggplot2::geom_point() +
        ggplot2::ylab("Predicted Z (95% interval)") +
        ggplot2::xlab("x along dimension i (other dims at random values)")
    },
    #' @description Plot leave one out predictions for design points
    # @importFrom ggplot2 ggplot aes stat_smooth geom_abline geom_segment
    # @importFrom ggplot2 geom_point geom_text xlab ylab ggtitle
    plotLOO = function() {
      ploo <- self$pred_LOO(se.fit = T)
      loodf <- cbind(ploo, Z=self$Z)
      loodf
      loodf$upper <- loodf$fit + 1.96 * loodf$se.fit
      loodf$lower <- loodf$fit - 1.96 * loodf$se.fit
      # Add text with coverage, R-sq
      coveragevec <- with(loodf, upper >= Z & lower <= Z)
      coverage <- mean(coveragevec)
      coverage
      rsq <- with(loodf, 1 - (sum((fit-Z)^2)) / (sum((mean(Z)-Z)^2)))
      rsq
      ggplot2::ggplot(loodf, ggplot2::aes(fit, Z)) +
        ggplot2::stat_smooth(method="loess", formula="y~x") +
        ggplot2::geom_abline(slope=1, intercept=0, color="red") +
        ggplot2::geom_segment(ggplot2::aes(x=lower, xend=upper, yend=Z), color="green") +
        ggplot2::geom_point() +
        # geom_text(x=min(loodf$fit), y=max(loodf$Z), label="abc") +
        ggplot2::geom_text(x=-Inf, y=Inf,
                           label=paste("Coverage:", signif(coverage,5)),
                           hjust=0, vjust=1) +
        ggplot2::geom_text(x=-Inf, y=Inf,
                           label=paste("R-sq:        ", signif(rsq,5)),
                           hjust=0, vjust=2.2) +
        # geom_text(x=Inf, y=-Inf, label="def", hjust=1, vjust=0)
        ggplot2::xlab("Predicted values (fit)") +
        ggplot2::ylab("Actual values (Z)") +
        ggplot2::ggtitle("Calibration of leave-one-out (LOO) predictions")
    },
    #' @description If track_optim, this will plot the parameters
    #' in the order they were evaluated.
    #' @param minindex Minimum index to plot.
    plot_track_optim = function(minindex=NULL) {
      if (length(self$track_optim_inputs) < 0.5) {
        stop("Can't plot_track_optim if track_optim was FALSE")
      }
      toi <- as.data.frame(matrix(unlist(self$track_optim_inputs), byrow=T,
                                  ncol=length(self$track_optim_inputs[[1]])))
      toi$deviance <- self$track_optim_dev
      toi$index <- 1:nrow(toi)
      if (!missing(minindex) && !is.null(minindex)) {
        stopifnot(is.numeric(minindex) && length(minindex) == 1)
        toi <- toi[toi$index >= minindex, ]
      }
      toi2 <- tidyr::pivot_longer(toi, cols = 1:(ncol(toi)-1))
      ggplot2::ggplot(toi2, ggplot2::aes(index, value)) +
        ggplot2::geom_line() +
        ggplot2::facet_wrap(.~name, scales='free_y')
    },
    #' @description Calculate loglikelihood of parameters
    #' @param mu Mean parameters
    #' @param s2 Variance parameter
    loglikelihood = function(mu=self$mu_hatX, s2=self$s2_hat) {
      -.5 * (self$N*log(s2) + as.numeric(determinant(self$K,logarithm=TRUE)$modulus) + #log(det(self$K)) +
               t(self$Z - mu)%*%self$Kinv%*%(self$Z - mu)/s2)
    },
    #' @description Get optimization functions
    #' @param param_update Should parameters be updated?
    #' @param nug.update Should nugget be updated?
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
      #       self$deviance_grad(params=params[1:(l-1)], nuglog=params[l],
      #                          nug.update=nug.update)
      #     },
      #     fngr=function(params) {
      #       l <- length(params)
      #       list(
      #         fn=function(params) {
      #           self$deviance(params=params[1:(l-1)], nuglog=params[l])
      #         },
      #         gr=function(params) {self$deviance_grad(
      #                 params=params[1:(l-1)], nuglog=params[l],
      #                 nug.update=nug.update)
      #         }
      #       )
      #     }
      #   )
      # } else {
      #   list(
      #     fn=function(params) {self$deviance(params=params)},
      #     gr=function(params) {self$deviance_grad(params=params,
      #                          nug.update=nug.update)},
      #     fngr=function(params) {
      #       list(
      #         fn=function(params) {self$deviance(params=params)},
      #         gr=function(params) {self$deviance_grad(params=params,
      #                              nug.update=nug.update)}
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
          if (self$track_optim) {
            self$track_optim_inputs[[length(self$track_optim_inputs)+1]] <- params
          }
          tparams <- if (tl>0) {params[ti]} else {NULL}
          kparams <- if (kl>0) {params[ki]} else {NULL}
          nparams <- if (nl>0) {params[ni]} else {NULL}
          dev <- self$deviance(params=kparams, nuglog=nparams,
                               trend_params=tparams)
          if (self$track_optim) {
            self$track_optim_dev[[length(self$track_optim_dev)+1]] <- dev
          }
          dev
        },
        gr=function(params) {
          if (self$track_optim) {
            self$track_optim_inputs[[length(self$track_optim_inputs)+1]] <- params
          }
          tparams <- if (tl>0) {params[ti]} else {NULL}
          kparams <- if (kl>0) {params[ki]} else {NULL}
          nparams <- if (nl>0) {params[ni]} else {NULL}
          dev_grad <- self$deviance_grad(params=kparams, nuglog=nparams,
                                         trend_params=tparams, nug.update=nug.update)
          if (self$track_optim) { # Doesn't actually get value
            self$track_optim_dev[[length(self$track_optim_dev)+1]] <- NA
          }
          dev_grad
        },
        fngr=function(params) {
          if (self$track_optim) {
            self$track_optim_inputs[[length(self$track_optim_inputs)+1]] <- params
          }
          tparams <- if (tl>0) {params[ti]} else {NULL}
          kparams <- if (kl>0) {params[ki]} else {NULL}
          nparams <- if (nl>0) {params[ni]} else {NULL}
          dev_fngr <- self$deviance_fngr(params=kparams, nuglog=nparams,
                                         trend_params=tparams, nug.update=nug.update)
          if (self$track_optim) {
            self$track_optim_dev[[length(self$track_optim_dev)+1]] <- dev_fngr$fn
          }
          dev_fngr
        }
      )
    },
    #' @description Lower bounds of parameters for optimization
    #' @param nug.update Is the nugget being updated?
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
    #' @description Upper bounds of parameters for optimization
    #' @param nug.update Is the nugget being updated?
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
    #' @description Starting point for parameters for optimization
    #' @param jitter Should there be a jitter?
    #' @param nug.update Is nugget being updated?
    param_optim_start = function(nug.update, jitter) {
      # param_start <- self$kernel$param_optim_start(jitter=jitter)
      if (nug.update) {
        nug_start <- log(self$nug,10)
        if (jitter) {
          nug_start <- nug_start + rexp(1, 1)
          nug_start <- min(max(log(self$nug.min,10),
                               nug_start),
                           log(self$nug.max,10))
        }
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
    #' @description Starting point for parameters for optimization
    #' @param jitter Should there be a jitter?
    #' @param nug.update Is nugget being updated?
    param_optim_start0 = function(nug.update, jitter) {
      # param_start <- self$kernel$param_optim_start0(jitter=jitter)
      if (nug.update) {
        nug_start <- -4
        if (jitter) {nug_start <- nug_start + rexp(1, 1)}
        # Make sure nug_start is in nug range
        nug_start <- min(max(log(self$nug.min,10), nug_start),
                         log(self$nug.max,10))
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
    #' @description Get matrix for starting points of optimization
    #' @param restarts Number of restarts to use
    #' @param nug.update Is nugget being updated?
    #' @param l Not used
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
    #' @description Optimize parameters
    #' @param restarts Number of restarts to do
    #' @param n0 This many starting parameters are chosen and evaluated.
    #' The best ones are used as the starting points for optimization.
    #' @param param_update Should parameters be updated?
    #' @param nug.update Should nugget be updated?
    #' @param parallel Should restarts be done in parallel?
    #' @param parallel_cores If running parallel, how many cores should be used?
    optim = function (restarts = self$restarts, n0=5*self$D,
                      param_update = T,
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
      n0 <- max(n0, restarts+1)
      param_optim_start_mat <- self$param_optim_start_mat(restarts=n0-1, #restarts,
                                                          nug.update=nug.update,
                                                          l=length(lower))
      if (!is.matrix(param_optim_start_mat)) {
        # Is a vector, should be a matrix with one row since it applies
        #    over columns
        param_optim_start_mat <- matrix(param_optim_start_mat, nrow=1)
      }
      # Below could just be run when condition is true,
      #  but it needs devs below anayways.
      if (TRUE || n0 > restarts + 1.5) {
        # Find best starting points
        devs <- rep(NA, ncol(param_optim_start_mat))
        for (i in 1:ncol(param_optim_start_mat)) {
          try(devs[i] <- optim.func(param_optim_start_mat[,i]), silent=T)
        }
        # Find best to start with
        best_start_inds <- order(order(devs))
        param_optim_start_mat <- param_optim_start_mat[, best_start_inds < restarts+1.5, drop=F]
      }


      # # This will make sure it at least can start
      # # Run before it sets initial parameters
      # # try.devlog <- try(devlog <- optim.func(start.par), silent = T)
      # try.devlog <- try(devlog <- optim.func(param_optim_start_mat[,1]),
      #                   silent = T)
      # if (inherits(try.devlog, "try-error") || is.infinite(devlog)) {
      #   warning("Current nugget doesn't work, increasing it #31973")
      #   # This will increase the nugget until cholesky works
      #   self$update_K_and_estimates()
      #   # devlog <- optim.func(start.par)
      #   devlog <- optim.func(param_optim_start_mat[,1])
      # }
      devlog <- devs[best_start_inds[which.min(best_start_inds)]]
      if (is.na(devlog) ||
          is.nan(devlog)) {
        warning("Current nugget doesn't work, increasing it #752983")
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
        message("All restarts had error, keeping initial")
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
    #' @param start.par.i Starting parameters for this restart
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
        {
          if (self$useGrad) {
            if (is.null(optim.fngr)) {
              lbfgs::lbfgs(optim.func, optim.grad, start.par.i, invisible=1)
            } else {
              # Two options for shared grad
              if (self$optimizer == "L-BFGS-B") {
                # optim uses L-BFGS-B which uses upper and lower
                # optim_share(fngr=optim.fngr, par=start.par.i,
                #             method='L-BFGS-B', upper=upper, lower=lower)
                # if (use_optim_share2) {
                optim_share2(fngr=optim.fngr, par=start.par.i,
                             method='L-BFGS-B', upper=upper, lower=lower)
                # } else {
                #   optim_share(fngr=optim.fngr, par=start.par.i,
                #               method='L-BFGS-B', upper=upper, lower=lower)
                # }
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
        }, silent=TRUE
      )
      if (!inherits(current, "try-error")) {
        # if (self$useGrad) {
        if (is.null(current$counts)) {current$counts <- c(NA,NA)}
        if(is.null(current$message)) {current$message=NA}
        # }
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
    },
    #' @description Fast update when adding new data.
    #' @param Xnew New X values to add.
    #' @param Znew New Z values to add.
    update_fast = function (Xnew=NULL, Znew=NULL) {
      # Updates data, K, and Kinv, quickly without adjusting parameters
      # Should be O(n^2) instead of O(n^3), but in practice not much faster
      N1 <- nrow(self$X)
      N2 <- nrow(Xnew)
      inds2 <- (N1+1):(N1+N2) # indices for new col/row, shorter than inds1
      K2 <- self$kernel$k(Xnew) + diag(self$kernel$s2 * self$nug, N2)
      K12 <- self$kernel$k(self$X, Xnew) # Need this before update_data

      self$update_data(Xnew=Xnew, Znew=Znew) # Doesn't update Kinv, etc

      # Update K
      K3 <- matrix(0, nrow=self$N, ncol=self$N)
      K3[-inds2, -inds2] <- self$K
      K3[-inds2, inds2] <- K12
      K3[inds2, -inds2] <- t(K12)
      K3[inds2, inds2] <- K2

      # Check for accuracy
      # summary(c(K3 - (self$kernel$k(self$X) +
      #                   diag(self$kernel$s2*self$nug, self$N))))

      self$K <- K3

      # Update the inverse using the block inverse formula
      K1inv_K12 <- self$Kinv %*% K12
      G <- solve(K2 - t(K12) %*% K1inv_K12)
      F1 <- -K1inv_K12 %*% G
      E <- self$Kinv - K1inv_K12 %*% t(F1)

      K3inv <- matrix(0, nrow=self$N, ncol=self$N)
      K3inv[-inds2, -inds2] <- E
      K3inv[-inds2, inds2] <- F1
      K3inv[inds2, -inds2] <- t(F1)
      K3inv[inds2, inds2] <- G

      # Check for accuracy
      # summary(c(K3inv - solve(self$K)))

      # self$K <- K3
      self$Kinv <- K3inv

      # self$mu_hatX <- self$trend$Z(X=self$X)
      # Just rbind new values
      self$mu_hatX <- rbind(self$mu_hatX,self$trend$Z(X=Xnew))

      invisible(self)
    },
    #' @description Update the parameters.
    #' @param nug.update Is the nugget being updated?
    #' @param ... Passed to optim.
    update_params = function(..., nug.update) {
      # Find lengths of params to optimize for each part
      tl <- length(self$trend$param_optim_start())
      kl <- length(self$kernel$param_optim_start())
      nl <- as.integer(nug.update)
      ti <- if (tl>0) {1:tl} else {c()}
      ki <- if (kl>0) {tl + 1:kl} else {c()}
      ni <- if (nl>0) {tl+kl+1:nl} else {c()}

      # If no params to optim, just return
      if (tl+kl+nl == 0) {return()}

      # Run optimization
      optim_out <- self$optim(..., nug.update=nug.update)

      # Split output into parts
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
    #' @description Update the data. Should only give in
    #' (Xnew and Znew) or (Xall and Zall).
    #' @param Xnew New X values to add.
    #' @param Znew New Z values to add.
    #' @param Xall All X values to be used. Will replace existing X.
    #' @param Zall All Z values to be used. Will replace existing Z.
    update_data = function(Xnew=NULL, Znew=NULL, Xall=NULL, Zall=NULL) {
      if (!is.null(Xall)) {
        self$X <- if (is.matrix(Xall)) Xall else matrix(Xall,nrow=1)
        self$N <- nrow(self$X)
      } else if (!is.null(Xnew)) {
        self$X <- rbind(self$X,
                        if (is.matrix(Xnew)) Xnew else matrix(Xnew,nrow=1))
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
        if (self$normalize) {
          Znewmat <- (Znewmat - self$normalize_mean) / self$normalize_sd
        }
        self$Z <- rbind(self$Z, Znewmat)
      }
      #if (!is.null(Xall) | !is.null(Xnew)) {self$update_K_and_estimates()}
      # update Kinv, etc, DONT THINK I NEED IT
    },
    #' @description Update correlation parameters. Not the nugget.
    #' @param ... Passed to self$update()
    update_corrparams = function (...) {
      self$update(nug.update = F, ...=...)
    },
    #' @description Update nugget Not the correlation parameters.
    #' @param ... Passed to self$update()
    update_nugget = function (...) {
      self$update(param_update = F, ...=...)
    },
    # deviance_searchnug = function() {
    #   optim(self$nug, function(nnug) {self$deviance(nug=nnug)},
    #       method="L-BFGS-B", lower=0, upper=Inf, hessian=F)$par
    # },
    # nugget_update = function () {
    #   nug <- self$deviance_searchnug()
    #   self$nug <- nug
    #   self$update_K_and_estimates()
    # },
    #' @description Calculate the deviance.
    #' @param params Kernel parameters
    #' @param nug Nugget
    #' @param nuglog Log of nugget. Only give in nug or nuglog.
    #' @param trend_params Parameters for the trend.
    deviance = function(params=NULL, nug=self$nug, nuglog, trend_params=NULL) {
      if (!missing(nuglog) && !is.null(nuglog)) {
        nug <- 10^nuglog
      }
      if (any(is.nan(params), is.nan(nug))) {
        if (self$verbose >= 2) {
          print("In deviance, returning Inf #92387")
        }
        return(Inf)
      }
      K <- self$kernel$k(x=self$X, params=params) +
        diag(nug, self$N) * self$kernel$s2_from_params(params=params)
      # if (is.nan(log(det(K)))) {return(Inf)}
      Z_hat <- self$trend$Z(X=self$X, params=trend_params)
      # dev.try <- try(dev <- log(det(K)) + sum((self$Z - self$mu_hat) *
      #                            solve(K, self$Z - self$mu_hat)))
      dev.try <- try(
        # dev <- log(det(K)) + sum((self$Z - Z_hat) * solve(K, self$Z - Z_hat))
        # Less likely to overflow by telling it to do logarithm
        dev <- (as.numeric(determinant(K,logarithm=TRUE)$modulus) +
                  sum((self$Z - Z_hat) * solve(K, self$Z - Z_hat)))

      )
      if (inherits(dev.try, "try-error")) {
        if (self$verbose>=2) {
          print("Deviance error #87126, returning Inf")
        }
        return(Inf)
      }
      # print(c(params, nuglog, dev))
      if (is.infinite(abs(dev))) {
        if (self$verbose>=2) {
          print("Deviance infinite #2332, returning Inf")
        }
        return(Inf)
      }
      dev
    },
    #' @description Calculate the gradient of the deviance.
    #' @param params Kernel parameters
    #' @param kernel_update Is the kernel being updated? If yes,
    #' it's part of the gradient.
    #' @param X Input matrix
    #' @param nug Nugget
    #' @param nug.update Is the nugget being updated? If yes,
    #' it's part of the gradient.
    #' @param nuglog Log of the nugget.
    #' @param trend_params Trend parameters
    #' @param trend_update Is the trend being updated? If yes,
    #' it's part of the gradient.
    deviance_grad = function(params=NULL, kernel_update=TRUE,
                             X=self$X,
                             nug=self$nug, nug.update, nuglog,
                             trend_params=NULL, trend_update=TRUE) {
      if (!missing(nuglog) && !is.null(nuglog)) {
        nug <- 10^nuglog
      }
      if (any(is.nan(params), is.nan(nug))) {
        if (self$verbose>=2) {
          print("In deviance_grad, returning NaN #92387")
        };
        return(rep(NaN, length(params)+as.integer(isTRUE(nug.update))))
      }
      C_nonug <- self$kernel$k(x=X, params=params)
      s2_from_kernel <- self$kernel$s2_from_params(params=params)
      C <- C_nonug + s2_from_kernel * diag(nug, self$N)
      # if (length(params) <= 0.5) {
      #   # Avoid error when no params are being updated
      #   kernel_update <- FALSE
      # } else {
      dC_dparams_out <- self$kernel$dC_dparams(params=params, X=X, C=C,
                                               C_nonug=C_nonug, nug=nug)
      dC_dparams <- dC_dparams_out#[[1]]
      # }
      # First of list should be list of dC_dparams
      # s2_from_kernel <- dC_dparams_out[[2]]
      # Second should be s2 for nugget deriv
      Z_hat <- self$trend$Z(X=X, params=trend_params)
      dZ_dparams <- self$trend$dZ_dparams(X=X, params=trend_params)
      # yminusmu <- self$Z - self$mu_hat
      yminusmu <- self$Z - Z_hat
      solve.try <- try(Cinv_yminusmu <- solve(C, yminusmu))
      if (inherits(solve.try, "try-error")) {
        if (self$verbose>=2) {
          print("Deviance grad error #63466, returning Inf")
        }
        return(Inf)
      }


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
    #' @description Calculate the deviance along with its gradient.
    #' @param params Kernel parameters
    #' @param kernel_update Is the kernel being updated? If yes,
    #' it's part of the gradient.
    #' @param X Input matrix
    #' @param nug Nugget
    #' @param nug.update Is the nugget being updated? If yes,
    #' it's part of the gradient.
    #' @param nuglog Log of the nugget.
    #' @param trend_params Trend parameters
    #' @param trend_update Is the trend being updated? If yes,
    #' it's part of the gradient.
    deviance_fngr = function(params=NULL, kernel_update=TRUE,
                             X=self$X,
                             nug=self$nug, nug.update, nuglog,
                             trend_params=NULL, trend_update=TRUE) {
      if (self$verbose >= 20) {cat('in deviance_fngr', '\n')}
      if (!missing(nuglog) && !is.null(nuglog)) {
        nug <- 10^nuglog
      }
      if (any(is.nan(params), is.nan(nug))) {
        if (self$verbose>=2) {
          print("In deviance_grad, returning NaN #92387")
        };
        return(rep(NaN, length(params)+as.integer(isTRUE(nug.update))))
      }
      # C_nonug <- self$kernel$k(x=X, params=params)
      # C <- C_nonug + s2_from_kernel * diag(nug, self$N)

      # s2_from_kernel <- self$kernel$s2_from_params(params=params)
      C_dC_try <- try(
        C_dC_dparams_out <- self$kernel$C_dC_dparams(params=params,
                                                     X=X, nug=nug),
        #C=C, C_nonug=C_nonug)
        silent = TRUE
      )
      if (inherits(C_dC_try, 'try-error')) {
        return(list(fn=self$deviance(params=params, nug=nug),
                    gr=self$deviance_grad(params=params, X=X, nug=nug,
                                          nug.update=nug.update)))
      }
      if (length(C_dC_dparams_out) < 2) {stop("Error #532987")}
      C <- C_dC_dparams_out[[1]]
      # First of list should be list of dC_dparams
      dC_dparams <- C_dC_dparams_out[[2]]
      # Second should be s2 for nugget deriv
      # s2_from_kernel <- dC_dparams_out[[2]]
      Z_hat <- self$trend$Z(X=X, params=trend_params)
      dZ_dparams <- self$trend$dZ_dparams(X=X, params=trend_params)
      # yminusmu <- self$Z - self$mu_hat
      yminusmu <- self$Z - Z_hat
      s2_from_kernel <- self$kernel$s2_from_params(params=params)
      # I tried not doing chol2inv, but it's faster inside of gradfunc, so keep it
      # if (calc_inv_from_chol) {
      Cinv <- chol2inv(chol(C))
      # solve.try <- try(Cinv_yminusmu <- solve(C, yminusmu))
      solve.try <- try(Cinv_yminusmu <- Cinv %*% yminusmu)
      # } else {
      #   chol_C <- chol(C)
      #   solve.try <- try(Cinv_yminusmu <- solvewithchol(chol_C, yminusmu))
      # }
      if (inherits(solve.try, "try-error")) {
        if (self$verbose>=2) {
          print("Deviance grad error #63466, returning Inf")
        }
        return(Inf)
      }

      gr <- c()
      if (length(dZ_dparams) > 0 && trend_update) {
        trend_gradfunc <- function(di) {
          # -2 * t(yminusmu) %*% solve(C, di) # Siginv %*% du/db
          # if (calc_inv_from_chol) {
          -2 * t(yminusmu) %*% (Cinv %*% di) # Siginv %*% du/db
          # } else {
          #   -2 * t(yminusmu) %*% solvewithchol(chol_C, di)
          # }
        }
        trend_gr <- apply(dZ_dparams, 2, trend_gradfunc)
        gr <- trend_gr
      } else {
        # trend_out <- c()
      }

      gradfunc <- function(di) {
        # t1 <- sum(diag(solve(C, di))) # Waste to keep resolving
        # t1 <- sum(diag((Cinv %*% di))) # Don't need whole mat mul
        # This is the main reason to calculate Cinv. Getting this trace this way
        # is way faster than having to repeatedly do solves with chol_C and then
        # taking the trace.
        t1 <- sum(Cinv * t(di))
        t2 <- sum(Cinv_yminusmu * (di %*% Cinv_yminusmu))
        t1 - t2
      }
      # out <- c(sapply(dC_dparams[[1]],gradfunc), gradfunc(dC_dparams[[2]]))
      # gr <- sapply(dC_dparams,gradfunc)
      if (kernel_update) {
        # Using apply() is 5x faster than Cpp code I wrote to do same thing
        #  Speed up by saving Cinv above to reduce number of solves
        # kernel_gr <- apply(dC_dparams, 1, gradfunc) # 6x faster below
        if (self$useC) {
          kernel_gr <- gradfuncarray(dC_dparams, Cinv, Cinv_yminusmu)
        } else {
          # Changing to R code so it works on my laptop
          kernel_gr <- gradfuncarrayR(dC_dparams, Cinv, Cinv_yminusmu)
        }
        gr <- c(gr, kernel_gr)
      }
      if (nug.update) {
        gr <- c(gr, gradfunc(diag(s2_from_kernel*nug*log(10), nrow(C))))
        # out <- c(out, gradfunc(diag(s2_from_kernel*, nrow(C)))*nug*log(10))
      }

      # Calculate fn
      logdetC <- as.numeric(determinant(C,logarithm=TRUE)$modulus) #log(det(C))
      if (is.nan(logdetC)) {
        dev <- Inf #return(Inf)
      } else {
        # dev.try <- try(dev <- logdetC + sum((yminusmu) * solve(C, yminusmu)))
        dev.try <- try(dev <- logdetC + sum((yminusmu) * Cinv_yminusmu))
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
            message(paste0("Deviance infinite #2333, returning 1e100, ",
                           "this is a hack and gives noticeable worse ",
                           "results on this restart."))
          }
          dev <- 1e100 # .Machine$double.xmax # Inf
        }
      }
      # dev

      # print(c(params, nuglog, out))
      out <- list(fn=dev, gr=gr)
      # cat('finished fngr, dev=', dev, ' par was', params, '\n')
      out
    },
    #' @description Calculate gradient
    #' @param XX points to calculate at
    #' @param X X points
    #' @param Z output points
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
      # Cinv_Z_minus_Zhat <- solve(self$K, Z - trendX)
      # Speed up since already have Kinv
      Cinv_Z_minus_Zhat <- self$Kinv %*% (Z - trendX)

      # Faster multiplication with arma_mult_cube_vec by 10x for big
      #  and .5x for small
      # t2 <- apply(dC_dx, 1, function(U) {U %*% Cinv_Z_minus_Zhat})
      t2 <- arma_mult_cube_vec(dC_dx, Cinv_Z_minus_Zhat)

      # if (ncol(dtrend_dx) > 1) {
      dtrend_dx + t(t2)
      # } else { # 1D needed transpose, not anymore
      # dtrend_dx + t2 # No longer needed with arma_mult_cube_vec
      # }
      # dtrend_dx + dC_dx %*% solve(self$K, Z - trendX)
    },
    #' @description Calculate norm of gradient
    #' @param XX points to calculate at
    grad_norm = function (XX) {
      grad1 <- self$grad(XX)
      if (!is.matrix(grad1)) return(abs(grad1))
      apply(grad1,1, function(xx) {sqrt(sum(xx^2))})
    },
    #' @description Calculate distribution of gradient
    #' @param XX points to calculate at
    grad_dist = function(XX) {
      # Calculates distribution of gradient at rows of XX
      if (!is.matrix(XX)) {
        if (self$D == 1) XX <- matrix(XX, ncol=1)
        else if (length(XX) == self$D) XX <- matrix(XX, nrow=1)
        else stop('grad_dist input should be matrix')
      } else {
        if (ncol(XX) != self$D) {stop("Wrong dimension input")}
      }
      nn <- nrow(XX)
      # Get mean from self$grad
      mn <- self$grad(XX=XX)
      # Calculate covariance here
      cv <- array(data = NA_real_, dim = c(nn, self$D, self$D))
      # New way calculates c1 and c2 outside loop
      c2 <- self$kernel$d2C_dudv_ueqvrows(XX=XX)
      c1 <- self$kernel$dC_dx(XX=XX, X=self$X)
      for (i in 1:nn) {
        # c2 <- self$kernel$d2C_dudv(XX=XX[i,,drop=F], X=XX[i,,drop=F])
        # c1 <- self$kernel$dC_dx(XX=XX[i,,drop=F], X=self$X)
        tc1i <- c1[i,,] # 1D gives problem, only need transpose if D>1
        if (!is.null(dim(tc1i))) {tc1i <- t(tc1i)}
        cv[i, , ] <- c2[i,,] - c1[i,,] %*% (self$Kinv %*% tc1i)
      }
      list(mean=mn, cov=cv)
    },
    #' @description Sample gradient at points
    #' @param XX points to calculate at
    #' @param n Number of samples
    grad_sample = function(XX, n) {
      if (!is.matrix(XX)) {
        if (length(XX) == self$D) {XX <- matrix(XX, nrow=1)}
        else {stop("Wrong dimensions #12574")}
      }
      # if (nrow(XX) > 1) {return(apply(XX, 1, self$grad_sample))}
      if (nrow(XX) > 1) {stop("Only can do 1 grad sample at a time")}
      grad_dist <- self$grad_dist(XX=XX)
      grad_samp <- MASS::mvrnorm(n=n, mu = grad_dist$mean[1,],
                                 Sigma = grad_dist$cov[1,,])
      grad_samp
      # gs2 <- apply(gs, 1, . %>% sum((.)^2))
      # c(mean(1/gs2), var(1/gs2))
    },
    #' @description Calculate mean of gradient norm squared
    #' @param XX points to calculate at
    grad_norm2_mean = function(XX) {
      # Calculate mean of squared norm of gradient
      # XX is matrix of points to calculate it at
      # Twice as fast as use self$grad_norm2_dist(XX)$mean
      grad_dist <- self$grad_dist(XX=XX)
      sum_grad_dist_mean_2 <- rowSums(grad_dist$mean^2)
      # Use sapply to get trace of cov matrix, return sum
      sum_grad_dist_mean_2 + sapply(1:nrow(XX), function(i) {
        if (ncol(XX)==1 ) {
          grad_dist$cov[i,,]
        } else {
          sum(diag(grad_dist$cov[i,,]))
        }
      })
    },
    #' @description Calculate distribution of gradient norm squared
    #' @param XX points to calculate at
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
        # Don't need to invert it, just solve with SigmaRoot
        # SigmaInv_i <- solve(Sigma_i)
        # Using my own sqrt function since it is faster.
        # # SigmaInvRoot_i <- expm::sqrtm(SigmaInv_i)
        # SigmaInvRoot_i <- sqrt_matrix(mat=SigmaInv_i, symmetric = TRUE)
        SigmaRoot_i <- sqrt_matrix(mat=Sigma_i, symmetric=TRUE)
        eigen_i <- eigen(Sigma_i)
        P_i <- t(eigen_i$vectors)
        lambda_i <- eigen_i$values
        # testthat::expect_equal(t(P) %*% diag(eth$values) %*% (P), Sigma)
        #               # Should be equal
        # b_i <- P_i %*% SigmaInvRoot_i %*% mean_i
        b_i <- P_i %*% solve(SigmaRoot_i, mean_i)
        g2mean_i <- sum(lambda_i * (b_i^2 + 1))
        g2var_i <- sum(lambda_i^2 * (4*b_i^2+2))
        means[i] <- g2mean_i
        vars[i] <- g2var_i
      }
      data.frame(mean=means, var=vars)
    },
    #' @description Get samples of squared norm of gradient
    #' @param XX points to sample at
    #' @param n Number of samples
    grad_norm2_sample = function(XX, n) {
      # Get samples of squared norm of gradient, check with grad_norm2_dist
      d <- ncol(XX)
      nn <- nrow(XX)
      out_sample <- matrix(NA_real_, nn, n)
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
        # testthat::expect_equal(t(P) %*% diag(eth$values) %*%
        #                       (P), Sigma) # Should be equal
        b_i <- c(P_i %*% SigmaInvRoot_i %*% mean_i)
        out_sample[i, ] <- replicate(n,
                                     sum(lambda_i * (rnorm(d) + b_i) ^ 2))
      }
      out_sample
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
    #' @description Calculate Hessian
    #' @param XX Points to calculate Hessian at
    #' @param as_array Should result be an array?
    hessian = function(XX, as_array=FALSE) {
      if (!is.matrix(XX)) {
        if (self$D == 1) XX <- matrix(XX, ncol=1)
        else if (length(XX) == self$D) XX <- matrix(XX, nrow=1)
        else stop('Predict input should be matrix')
      } else {
        if (ncol(XX) != self$D) {stop("Wrong dimension input")}
      }
      hess1 <- array(NA_real_, dim = c(nrow(XX), self$D, self$D))
      for (i in 1:nrow(XX)) { # 0 bc assume trend has zero hessian
        d2 <- self$kernel$d2C_dx2(XX=XX[i,,drop=F], X=self$X)
        for (j in 1:self$D) {
          hess1[i, j, ] <- 0 + d2[1,j,,] %*% self$Kinv %*%
            (self$Z - self$mu_hatX)
        }
      }
      if (nrow(XX) == 1 && !as_array) { # Return matrix if only one value
        hess1[1,,]
      } else {
        hess1
      }
    },
    #' @description Sample at rows of XX
    #' @param XX Input matrix
    #' @param n Number of samples
    sample = function(XX, n=1) {
      # Generates n samples at rows of XX
      px <- self$pred(XX, covmat = T)
      Sigma.try <- try(newy <- MASS::mvrnorm(n=n, mu=px$mean, Sigma=px$cov))
      if (inherits(Sigma.try, "try-error")) {
        message("Adding nugget to get sample")
        Sigma.try2 <- try(newy <- MASS::mvrnorm(n=n, mu=px$mean,
                                                Sigma=px$cov +
                                                  diag(self$nug, nrow(px$cov))))
        if (inherits(Sigma.try2, "try-error")) {
          stop("Can't do sample, can't factor Sigma")
        }
      }
      newy # Not transposing matrix since it gives var a problem
    },
    #' @description Calculate expected improvement
    #' @param x Vector to calculate EI of, or matrix for whose rows it should
    #' be calculated
    #' @param minimize Are you trying to minimize the output?
    #' @param eps Exploration parameter
    EI = function(x, minimize=FALSE, eps=0) {
      stopifnot(length(minimize)==1, is.logical(minimize))
      stopifnot(length(eps)==1, is.numeric(eps), eps >= 0)
      # if (minimize) {
      #   stop('can only max for EI, not min')
      # }
      if (is.matrix(x)) {
        stopifnot(ncol(x) == ncol(self$X))
      } else if (is.vector(x)) {
        stopifnot(length(x) == ncol(self$X))
      } else {
        stop(paste0("bad x in EI, class is: ", class(x)))
      }
      # stopifnot(is.vector(x), length(x) == ncol(self$X))
      # fxplus <- if (minimize) {min(self$Z)} else {max(self$Z)}
      # pred <- self$pred(x, se.fit=T)
      # Need to use prediction of mean
      xnew_meanpred <- self$pred(x, se.fit=T, mean_dist=T)
      Xold_meanpred <- self$pred(self$X, se.fit=T, mean_dist=T)
      # Use predicted mean at each point since it doesn't make sense not to
      # when there is noise. Or should fxplus be optimized over inputs?
      fxplus <- if (minimize) {min(Xold_meanpred$mean)} else {max(Xold_meanpred$mean)}
      if (minimize) {
        # Ztop <- fxplus - pred$mean - eps
        Ztop <- fxplus - xnew_meanpred$mean - eps
      } else {
        # Ztop <- pred$mean - fxplus - eps
        Ztop <- xnew_meanpred$mean - fxplus - eps
      }
      # Z <- Ztop / pred$se
      Z <- Ztop / xnew_meanpred$se
      # if (pred$se <= 0) {return(0)}
      # (Ztop) * pnorm(Z) + pred$se * dnorm(Z)
      # ifelse(pred$se <= 0, 0,
      #        (Ztop) * pnorm(Z) + pred$se * dnorm(Z))
      ifelse(xnew_meanpred$se <= 0, 0,
             (Ztop) * pnorm(Z) + xnew_meanpred$se * dnorm(Z))
    },
    #' @description Find the point that maximizes the expected improvement
    #' @param lower Lower bounds to search within
    #' @param upper Upper bounds to search within
    #' @param n0 Number of points to evaluate in initial stage
    #' @param minimize Are you trying to minimize the output?
    #' @param eps Exploration parameter
    maxEI = function(lower=apply(self$X, 2, min), upper=apply(self$X, 2, max),
                     n0=100, minimize=FALSE, eps=0,
                     discreteinputs=NULL) {
      stopifnot(all(lower < upper))
      stopifnot(length(n0)==1, is.numeric(n0), n0>=1)
      # Check if any kernels have factors
      if (!is.null(discreteinputs) || length(find_kernel_factor_dims(self$kernel)) > 0) {
        # Has at least one factor
        # return(self$maxEIwithfactors(lower=lower, upper=upper, n0=n0,
        #                              minimize=minimize, eps=eps))
        return(self$maxEIwithfactorsordiscrete(lower=lower, upper=upper, n0=n0,
                                               minimize=minimize, eps=eps,
                                               discreteinputs=discreteinputs))
      }
      # Random points to evaluate to find best starting point
      X0 <- lhs::randomLHS(n=n0, k=ncol(self$X))
      X0 <- sweep(X0, 2, upper-lower, "*")
      X0 <- sweep(X0, 2, lower, "+")
      # Also use point near current optimum
      bestZind <- if (minimize) {which.min(self$Z)} else {which.max(self$Z)}
      X0 <- rbind(X0,
                  pmax(pmin(self$X[bestZind, ] +
                              rnorm(ncol(self$X), 0, 1e-4),
                            upper),
                       lower))
      # Calculate EI at these points
      EI0 <- self$EI(x=X0, minimize=minimize, eps=eps)
      if (all(EI0 <= 0)) {
        warning("maxEI couldn't find any inputs with nonzero EI")
      }
      ind <- which.max(EI0)
      # Optimize starting from that point to find input that maximizes EI
      optim_out <- optim(par=X0[ind,],
                         lower=lower, upper=upper,
                         # fn=function(xx){ei <- -self$EI(xx); cat(xx, ei, "\n"); ei},
                         fn=function(xx){-self$EI(xx, minimize = minimize)},
                         method="L-BFGS-B")
      optim_out$par
      # Return list, same format as DiceOptim::max_EI
      list(
        par=optim_out$par,
        value=-optim_out$val
      )
    },
    #' @description Find the point that maximizes the expected improvement.
    #' Used whenever one of the inputs is a factor (can only take values 1:n).
    #' @param lower Lower bounds to search within
    #' @param upper Upper bounds to search within
    #' @param n0 Number of points to evaluate in initial stage
    #' @param minimize Are you trying to minimize the output?
    #' @param eps Exploration parameter
    maxEIwithfactors = function(lower=apply(self$X, 2, min),
                                upper=apply(self$X, 2, max),
                                n0=100, minimize=FALSE, eps=0) {
      stopifnot(all(lower < upper))
      stopifnot(length(n0)==1, is.numeric(n0), n0>=1)
      # Get factor info
      factorinfo <- find_kernel_factor_dims(self$kernel)
      # Run inner EI over all factor combinations
      stopifnot(length(factorinfo)>0)
      # factordf <- data.frame(index=factorinfo[1]
      factorlist <- list()
      for (i_f in 1:(length(factorinfo)/2)) {
        factorlist[[as.character(factorinfo[i_f*2-1])]] <- 1:factorinfo[i_f*2]
      }
      factordf <- do.call(expand.grid, factorlist)
      # Track best seen in optimizing EI
      bestval <- Inf
      bestpar <- c()
      factorxindex <- factorinfo[(1:(length(factorinfo)/2))*2-1] #factorinfo[[1]]
      factornlevels <- factorinfo[(1:(length(factorinfo)/2))*2]

      # If no non-factor levels, just predict EI at all and return best
      if (length(factorxindex) == self$D) {
        Xmat <- as.matrix(factordf)
        EI_Xmat <- self$EI(Xmat, minimize = minimize)
        bestind <- which.max(EI_Xmat)[1]
        bestval <- EI_Xmat[bestind]
        bestpar <- Xmat[bestind, ]
        return(unname(bestpar))
      } else {
        # Has continuous and factor indices
        # Alternate optimizing over cts and factors
        X0 <- lhs::randomLHS(n=n0, k=self$D)
        X0 <- sweep(X0, 2, upper-lower, "*")
        X0 <- sweep(X0, 2, lower, "+")
        for (j in 1:length(factorxindex)) {
          X0[, factorxindex[j]] <- sample(1:factornlevels[j], n0, replace=TRUE)
        }
        # Calculate EI at these points, use best as starting point for optim
        EI0 <- self$EI(x=X0, minimize=minimize, eps=eps)
        ind <- which.max(EI0)

        # Continuous indexes
        ctsinds <- setdiff(1:self$D, factorxindex)
        Xstart <- X0[ind, ]
        Xstartfactors <- Xstart[factorxindex]
        bestEIsofar <- EI0[ind]
        notdone <- TRUE
        i_while <- 0
        while(notdone) {
          # cat('in while loop', i_while, bestEIsofar, "\n")
          i_while <- i_while + 1
          # Optimize over cts variables
          # Optimize starting from that point to find input that maximizes EI
          optim_out_i_indcomb <- optim(par=Xstart[-factorxindex], #X0[ind, -factorxindex],
                                       lower=lower[-factorxindex],
                                       upper=upper[-factorxindex],
                                       # fn=function(xx){ei <- -self$EI(xx); cat(xx, ei, "\n"); ei},
                                       fn=function(xx){
                                         xx2 <- numeric(self$D)
                                         xx2[ctsinds] <- xx
                                         xx2[factorxindex] <- Xstartfactors
                                         # xx2 <- c(xx[xxinds1], factorxlevel, xx[xxinds2])
                                         # cat(xx, xx2, "\n")
                                         -self$EI(xx2, minimize = minimize)
                                       },
                                       method="L-BFGS-B")
          Xstart2 <- Xstart
          Xstart2[-factorxindex] <- optim_out_i_indcomb$par
          # Optimize over factors
          Xmat <- matrix(Xstart2, byrow=TRUE, nrow=nrow(factordf), ncol=self$D)
          Xmat[, factorxindex] <- as.matrix(factordf)
          EI_Xmat <- self$EI(Xmat, minimize = minimize)
          # for (i in 1:nrow(Xmat)) {
          #   Xmat[i, -factorxindex] <-
          # }
          newbestEI <- max(EI_Xmat)
          newbestX <- Xmat[which.min(EI_Xmat)[1],]
          # If no improvement, end it
          if (newbestEI - bestEIsofar <=  1e-16) {
            # return(newbestX)
            # Return list, same format as DiceOptim::max_EI
            return(
              list(
                par=newbestX,
                value=newbestEI
              )
            )
          }
          # Or break if enough iterations
          if (i_while > 10) {
            # return(newbestX)
            # Return list, same format as DiceOptim::max_EI
            return(
              list(
                par=newbestX,
                value=newbestEI
              )
            )
          }
          bestEIsofar <- newbestEI
        }

      }
      # stopifnot(length(bestpar) == self$D)
      # return(bestpar)
      stop("maxEIwithfactors failed #320198471")
    },
    maxEIwithfactorsordiscrete = function(lower=apply(self$X, 2, min),
                                          upper=apply(self$X, 2, max),
                                          n0=100, minimize=FALSE, eps=0,
                                          discreteinputs=NULL) {
      stopifnot(all(lower < upper))
      stopifnot(length(n0)==1, is.numeric(n0), n0>=1)
      # Make sure discreteinputs is okay
      if (!is.null(discreteinputs)) {
        stopifnot(is.list(discreteinputs), length(discreteinputs)>0)
        stopifnot(!is.null(names(discreteinputs)),
                  all(names(discreteinputs) != ""))
        discreteinds <- as.integer(names(discreteinputs))
        stopifnot(!is.na(discreteinds))
      } else {
        discreteinds <- c()
      }
      # Get factor info
      factorinfo <- find_kernel_factor_dims(self$kernel)
      if (!is.null(factorinfo)) {
        # Run inner EI over all factor combinations
        stopifnot(length(factorinfo)>0)
        # factordf <- data.frame(index=factorinfo[1]
        factorlist <- list()
        for (i_f in 1:(length(factorinfo)/2)) {
          factorlist[[as.character(factorinfo[i_f*2-1])]] <- 1:factorinfo[i_f*2]
        }
        factordf <- do.call(expand.grid, factorlist)
        # Track best seen in optimizing EI
        bestval <- Inf
        bestpar <- c()
        factorxindex <- factorinfo[(1:(length(factorinfo)/2))*2-1] #factorinfo[[1]]
        factornlevels <- factorinfo[(1:(length(factorinfo)/2))*2]
      } else {
        # Indices of factor columns
        factorxindex <- c()
      }

      ctsinds <- setdiff(1:self$D, c(discreteinds, factorxindex))

      # browser("do coordinate descent here?")

      # If no non-factor levels, just predict EI at all and return best
      if (length(factorxindex) == self$D) {
        # browser()
        Xmat <- as.matrix(factordf)
        EI_Xmat <- self$EI(Xmat, minimize = minimize)
        bestind <- which.max(EI_Xmat)[1]
        bestval <- EI_Xmat[bestind]
        bestpar <- Xmat[bestind, ]
        return(unname(bestpar))
      } else {
        # Has continuous/factor/discrete indices
        # Alternate optimizing over cts and factors an discrete
        # browser()
        #
        X0 <- lhs::randomLHS(n=n0, k=self$D)
        X0 <- sweep(X0, 2, upper-lower, "*")
        X0 <- sweep(X0, 2, lower, "+")
        for (j in seq_along(factorxindex)) {
          X0[, factorxindex[j]] <- sample(1:factornlevels[j], n0, replace=TRUE)
        }
        for (j in seq_along(discreteinds)) {
          X0[, discreteinds[j]] <- sample(discreteinputs[[j]], n0, replace=TRUE)
        }
        # Calculate EI at these points, use best as starting point for optim
        EI0 <- self$EI(x=X0, minimize=minimize, eps=eps)
        ind <- which.max(EI0)

        # Continuous indexes
        # ctsinds <- setdiff(1:self$D, factorxindex)
        Xstart <- X0[ind, ]
        # Xstartfactors <- Xstart[factorxindex]
        bestEIsofar <- EI0[ind]
        notdone <- TRUE
        # browser()
        i_while <- 0
        while(notdone) {
          # cat('in while loop', i_while, bestEIsofar, "\n")
          i_while <- i_while + 1
          # Optimize over cts variables
          if (length(ctsinds) > 0) {
            # Optimize starting from that point to find input that maximizes EI
            optim_out_i_indcomb <- optim(par=Xstart[ctsinds], #X0[ind, -factorxindex],
                                         lower=lower[ctsinds],
                                         upper=upper[ctsinds],
                                         # fn=function(xx){ei <- -self$EI(xx); cat(xx, ei, "\n"); ei},
                                         fn=function(xx){
                                           xx2 <- numeric(self$D)
                                           xx2[ctsinds] <- xx
                                           xx2[-ctsinds] <- Xstart[-ctsinds]
                                           # xx2 <- c(xx[xxinds1], factorxlevel, xx[xxinds2])
                                           # cat(xx, xx2, "\n")
                                           -self$EI(xx2, minimize = minimize)
                                         },
                                         method="L-BFGS-B")
            # Xstart2 <- Xstart
            # Xstart2[-factorxindex] <- optim_out_i_indcomb$par
            Xstart[ctsinds] <- optim_out_i_indcomb$par
            newbestEI <- -optim_out_i_indcomb$value
          }
          # Optimize over factors
          if (length(factorxindex) > 0) {
            Xmat <- matrix(Xstart, byrow=TRUE, nrow=nrow(factordf), ncol=self$D)
            Xmat[, factorxindex] <- as.matrix(factordf)
            EI_Xmat <- self$EI(Xmat, minimize = minimize)
            # for (i in 1:nrow(Xmat)) {
            #   Xmat[i, -factorxindex] <-
            # }
            bestfactorind <- which.max(EI_Xmat)[1]
            Xstart[factorxindex] <- Xmat[bestfactorind, factorxindex]
            newbestEI <- max(EI_Xmat)
          }
          # Optimize over discrete
          if (length(discreteinds) > 0) {
            for (i in seq_along(discreteinds)) {
              # ndiscrete <- length(discreteinputs)
              discretevals_i <- discreteinputs[[i]]
              # If too many, sample plus keep current best
              if (length(discretevals_i) > 1100) {
                discretevals_i <- c(Xstart[discreteinds[i]],
                                    sample(discretevals_i, 1000, replace=FALSE))
              }
              Xmat <- matrix(Xstart, byrow=TRUE,
                             nrow=length(discretevals_i),
                             ncol=self$D)
              Xmat[, discreteinds[i]] <- discretevals_i
              EI_Xmat <- self$EI(Xmat, minimize = minimize)
              bestdiscreteind <- which.max(EI_Xmat)[1]
              Xstart[discreteinds] <- Xmat[bestdiscreteind, discreteinds]
              newbestEI <- max(EI_Xmat)
            }
          }

          # newbestEI <- max(EI_Xmat)
          # newbestX <- Xmat[which.min(EI_Xmat)[1],]
          # If no improvement, end it
          if (newbestEI - bestEIsofar <= 1e-16) {
            # return(Xstart)
            # Return list, same format as DiceOptim::max_EI
            return(
              list(
                par=Xstart,
                value=newbestEI
              )
            )
          }
          # browser()
          # Or break if enough iterations
          if (i_while > 10) {
            # return(Xstart)
            # Return list, same format as DiceOptim::max_EI
            return(
              list(
                par=Xstart,
                value=newbestEI
              )
            )
          }
          bestEIsofar <- newbestEI
        }

      }
      # stopifnot(length(bestpar) == self$D)
      # return(bestpar)
      stop("maxEIwithfactorsordiscrete failed #259287")
    },
    #' @description Find the multiple points that maximize the expected
    #' improvement. Currently only implements the constant liar method.
    #' @param npoints Number of points to add
    #' @param method Method to use. Can only be "CL" for constant liar.
    #' @param lower Lower bounds to search within
    #' @param upper Upper bounds to search within
    #' @param n0 Number of points to evaluate in initial stage
    #' @param minimize Are you trying to minimize the output?
    #' @param eps Exploration parameter
    maxqEI = function(npoints, method="CL",
                      lower=apply(self$X, 2, min), upper=apply(self$X, 2, max),
                      n0=100, minimize=FALSE, eps=0) {
      stopifnot(is.numeric(npoints), length(npoints)==1, npoints >= 1)
      if (npoints==1) {
        # For single point, use proper function
        return(self$maxEI(lower=lower, upper=upper, n0=n0,
                          minimize=minimize, eps=eps))
      }
      stopifnot(method %in% c("CL", "pred"))
      # Clone object since we will add fake data
      gpclone <- self$clone(deep=TRUE)
      # Track points selected
      selectedX <- matrix(data=NA, nrow=npoints, ncol=ncol(self$X))
      Xmeanpred <- self$pred(self$X, se.fit=T, mean_dist=T)
      # Zimpute <- if (minimize) {min(self$Z)} else {max(self$Z)}
      # Constant liar value
      Zimpute <- if (minimize) {min(Xmeanpred$mean)} else {max(Xmeanpred$mean)}
      for (i in 1:npoints) {
        # Find and store point that maximizes EI
        maxEI_i <- gpclone$maxEI(lower=lower, upper=upper, n0=n0, eps=eps, minimize=minimize)
        xi <- maxEI_i$par
        selectedX[i, ] <- xi
        if (method == "pred") {
          # browser()
          Zimpute <- self$predict(xi)
        }
        # Update clone with new data, don't update parameters since it's fake data
        if (i < npoints) {
          gpclone$update(Xnew=xi, Znew=Zimpute, no_update=TRUE)
        }
      }
      # Return matrix of points
      # selectedX

      # Return list, same format as DiceOptim::max_EI
      return(
        list(
          par=selectedX,
          value=NA
        )
      )
    },
    KG = function(x, minimize=FALSE, eps=0, current_extreme=NULL) {
      if (exists('kgbrow') && kgbrow) {browser()}
      xkg <- x
      if (!is.matrix(xkg)) {
        stopifnot(length(xkg) == self$D)
        xkg <- matrix(xkg, nrow=1)
      }
      if (missing(current_extreme) || is.null(current_extreme)) {
        # Find current max/min
        # Find current max
        gpkgmax <- optim(par=self$X[which.max(self$Z)[1],],
                         fn=function(xx) {-self$pred(xx)},
                         method='Brent', lower=0, upper=1)
        current_extreme <- -gpkgmax$value
      } else {
        stopifnot(is.numeric(current_extreme), length(current_extreme) == 1)
      }
      # Sample at xkg
      xkgpred <- gpkg$pred(xkg, se.fit = T)
      xkgpred
      nsamps <- 5
      xkgsamps <- qnorm(((1:nsamps)-.5)/nsamps, xkgpred$mean, xkgpred$se)
      kgs <- rep(NA, nsamps)
      gpkgclone <- gpkg$clone(deep=TRUE)
      for (i in 1:nsamps) {
        xkgsamp <- xkgsamps[i]
        # xkgsamp <- rnorm(1, xkgpred$mean, xkgpred$se)
        # Add samp to mod
        # gpkgclone <- gpkg$clone(deep=TRUE)
        # gpkgclone$update(Xnew=xkg, Znew=xkgsamp, no_update = TRUE)
        gpkgclone$update(Xall=rbind(self$X, xkg),
                         Zall=rbind(self$Z, xkgsamp),
                         no_update = TRUE)
        # gpkgclone$plot1D()
        # Find clone max after adding sample
        gpkgmaxclone <- optim(par=gpkgclone$X[which.max(gpkgclone$Z)[1],],
                              fn=function(xx) {-gpkgclone$pred(xx)},
                              method='Brent', lower=0, upper=1)
        gpkgmaxclone
        # gpkgmaxclone$value - gpkgmax$value
        kgs[i] <- (-gpkgmaxclone$value) - current_extreme #gpkgmax$value
      }
      kgs
      mean(kgs)
    },
    #' @description Print this object
    print = function() {
      cat("GauPro kernel model object\n")
      cat(paste0("\tD = ", self$D, ", N = ", self$N,"\n"))
      cat(paste0("\tNugget = ", signif(self$nug, 3), "\n"))
      cat("\tRun $update() to add data and/or optimize again\n")
      cat("\tUse $pred() to get predictions at new points\n")
      cat("\tUse $plot() to visualize the model\n")
      invisible(self)
    }
  ),
  private = list(

  )
)
