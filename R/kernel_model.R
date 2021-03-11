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
        #' @param ... Not used
        initialize = function(X, Z,
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
        #' @description Predict for a matrix of points
        #' @param XX points to predict at
        #' @param se.fit Should standard error be returned?
        #' @param covmat Should covariance matrix be returned?
        #' @param split_speed Should the matrix be split for faster predictions?
        predict = function(XX, se.fit=F, covmat=F, split_speed=F) {
          self$pred(XX=XX, se.fit=se.fit, covmat=covmat,
                    split_speed=split_speed)
        },
        #' @description Predict for a matrix of points
        #' @param XX points to predict at
        #' @param se.fit Should standard error be returned?
        #' @param covmat Should covariance matrix be returned?
        #' @param split_speed Should the matrix be split for faster predictions?
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
        #' @description Predict for a matrix of points
        #' @param XX points to predict at
        #' @param se.fit Should standard error be returned?
        #' @param covmat Should covariance matrix be returned?
        #' @param return_df When returning se.fit, should it be returned in
        #' a data frame?
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
          Sigma.try <- try(newy <- MASS::mvrnorm(n=n2, mu=px$mean,
                                                 Sigma=px$cov),
                           silent = TRUE)
          if (inherits(Sigma.try, "try-error")) {
            message("Adding nugget to cool1Dplot")
            Sigma.try2 <- try(
              newy <- MASS::mvrnorm(n=n2, mu=px$mean,
                                  Sigma=px$cov + diag(self$nug, nrow(px$cov))))
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
                 if (self$normalize) {
                   self$Z * self$normalize_sd + self$normalize_mean
                 } else {self$Z},
                 pch=19, col=1, cex=2)
        },
        #' @description Make 1D plot
        #' @param n2 Number of things to plot
        #' @param nn Number of things to plot
        #' @param col2 color
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
                                    pts=self$X)
        },
        #' @description Calculate loglikelihood of parameters
        #' @param mu Mean parameters
        #' @param s2 Variance parameter
        loglikelihood = function(mu=self$mu_hatX, s2=self$s2_hat) {
          -.5 * (self$N*log(s2) + log(det(self$K)) +
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
              tparams <- if (tl>0) {params[ti]} else {NULL}
              kparams <- if (kl>0) {params[ki]} else {NULL}
              nparams <- if (nl>0) {params[ni]} else {NULL}
              self$deviance(params=kparams, nuglog=nparams,
                            trend_params=tparams)
            },
            gr=function(params) {
              tparams <- if (tl>0) {params[ti]} else {NULL}
              kparams <- if (kl>0) {params[ki]} else {NULL}
              nparams <- if (nl>0) {params[ni]} else {NULL}
              self$deviance_grad(params=kparams, nuglog=nparams,
                                 trend_params=tparams, nug.update=nug.update)
            },
            fngr=function(params) {
              tparams <- if (tl>0) {params[ti]} else {NULL}
              kparams <- if (kl>0) {params[ki]} else {NULL}
              nparams <- if (nl>0) {params[ni]} else {NULL}
              self$deviance_fngr(params=kparams, nuglog=nparams,
                                 trend_params=tparams, nug.update=nug.update)
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
        #' @param param_update Should parameters be updated?
        #' @param nug.update Should nugget be updated?
        #' @param parallel Should restarts be done in parallel?
        #' @param parallel_cores If running parallel, how many cores should be used?
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
          if (is.nan(log(det(K)))) {browser();return(Inf)}
          Z_hat <- self$trend$Z(X=self$X, params=trend_params)
          # dev.try <- try(dev <- log(det(K)) + sum((self$Z - self$mu_hat) *
          #                            solve(K, self$Z - self$mu_hat)))
          dev.try <- try(
            dev <- log(det(K)) + sum((self$Z - Z_hat) * solve(K, self$Z - Z_hat))
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
          dC_dparams_out <- self$kernel$dC_dparams(params=params, X=X, C=C,
                                                   C_nonug=C_nonug, nug=nug)
          dC_dparams <- dC_dparams_out#[[1]]
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
                                 trend_params=NULL, trend_update=TRUE) {#browser()
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
          Cinv <- chol2inv(chol(C))
          # solve.try <- try(Cinv_yminusmu <- solve(C, yminusmu))
          solve.try <- try(Cinv_yminusmu <- Cinv %*% yminusmu)
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
              -2 * t(yminusmu) %*% (Cinv %*% di) # Siginv %*% du/db
            }
            trend_gr <- apply(dZ_dparams, 2, trend_gradfunc)
            gr <- trend_gr
          } else {
            # trend_out <- c()
          }

          gradfunc <- function(di) {
            # t1 <- sum(diag(solve(C, di))) # Waste to keep resolving
            # t1 <- sum(diag((Cinv %*% di))) # Don't need whole mat mul
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
            kernel_gr <- gradfuncarray(dC_dparams, Cinv, Cinv_yminusmu)
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
                print("Deviance infinite #2333, returning 1e100,
                      this is a hack and gives noticeable worse
                      results on this restart.")
              }
              dev <- 1e100 # .Machine$double.xmax # Inf
            }
          }
          # dev

          # print(c(params, nuglog, out))
          out <- list(fn=dev, gr=gr)
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
        hessian = function(XX, as_array=FALSE) {#browser()
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
        #' @description Print this object
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
