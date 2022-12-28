#' Latent Factor Kernel R6 class
#'
#' Used for factor variables, a single dimension.
#' Each level of the factor gets mapped into a latent space,
#' then the distances in that space determine their correlations.
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @useDynLib GauPro, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats optim
# @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @field p Parameter for correlation
#' @field p_est Should p be estimated?
#' @field p_lower Lower bound of p
#' @field p_upper Upper bound of p
#' @field p_length length of p
#' @field s2 variance
#' @field s2_est Is s2 estimated?
#' @field logs2 Log of s2
#' @field logs2_lower Lower bound of logs2
#' @field logs2_upper Upper bound of logs2
#' @field xindex Index of the factor (which column of X)
#' @field nlevels Number of levels for the factor
#' @field latentdim Dimension of embedding space
#' @field pf_to_p_log Logical vector used to convert pf to p
#' @field p_to_pf_inds Vector of indexes used to convert p to pf
#' @field offdiagequal What should offdiagonal values be set to when the
#' indices are the same? Use to avoid decomposition errors, similar to
#' adding a nugget.
#' @examples
#' # Create a new kernel for a single factor with 5 levels,
#' #  mapped into two latent dimensions.
#' kk <- LatentFactorKernel$new(D=1, nlevels=5, xindex=1, latentdim=2)
#' # Random initial parameter values
#' kk$p
#' # Plots to understand
#' kk$plotLatent()
#' kk$plot()
# kmat <- outer(1:5, 1:5, Vectorize(kk$k))
# kmat
#'
# kk$dC_dparams(X=matrix(1:5, ncol=1), nug=0)
# kk$C_dC_dparams(X=matrix(1:5, ncol=1), nug=0, params=c(kk$p, kk$s2))$C
#'
#' # 5 levels, 1/4 are similar and 2/3/5 are similar
#' n <- 30
#' x <- matrix(sample(1:5, n, TRUE))
#' y <- c(ifelse(x == 1 | x == 4, 4, -3) + rnorm(n,0,.1))
#' plot(c(x), y)
#' m5 <- GauPro_kernel_model$new(
#'   X=x, Z=y,
#'   kernel=LatentFactorKernel$new(D=1, nlevels = 5, xindex = 1, latentdim = 2))
#' m5$kernel$p
#' # We should see 1/4 and 2/3/4 in separate clusters
#' m5$kernel$plotLatent()
#'
# 2D, Gaussian on 1D, LatentFactor on 2nd dim
#' library(dplyr)
#' n <- 20
#' X <- cbind(matrix(runif(n,2,6), ncol=1),
#'            matrix(sample(1:2, size=n, replace=TRUE), ncol=1))
#' X <- rbind(X, c(3.3,3), c(3.7,3))
#' n <- nrow(X)
#' Z <- X[,1] - (4-X[,2])^2 + rnorm(n,0,.1)
#' plot(X[,1], Z, col=X[,2])
#' tibble(X=X, Z) %>% arrange(X,Z)
#' k2a <- IgnoreIndsKernel$new(k=Gaussian$new(D=1), ignoreinds = 2)
#' k2b <- LatentFactorKernel$new(D=2, nlevels=3, xind=2, latentdim=2)
#' k2 <- k2a * k2b
#' k2b$p_upper <- .65*k2b$p_upper
#' gp <- GauPro_kernel_model$new(X=X, Z=Z, kernel = k2, verbose = 5,
#'   nug.min=1e-2, restarts=1)
#' gp$kernel$k1$kernel$beta
#' gp$kernel$k2$p
#' gp$kernel$k(x = gp$X)
#' tibble(X=X, Z=Z, pred=gp$predict(X)) %>% arrange(X, Z)
#' tibble(X=X[,2], Z) %>% group_by(X) %>% summarize(n=n(), mean(Z))
#' curve(gp$pred(cbind(matrix(x,ncol=1),1)),2,6, ylim=c(min(Z), max(Z)))
#' points(X[X[,2]==1,1], Z[X[,2]==1])
#' curve(gp$pred(cbind(matrix(x,ncol=1),2)), add=TRUE, col=2)
#' points(X[X[,2]==2,1], Z[X[,2]==2], col=2)
#' curve(gp$pred(cbind(matrix(x,ncol=1),3)), add=TRUE, col=3)
#' points(X[X[,2]==3,1], Z[X[,2]==3], col=3)
#' legend(legend=1:3, fill=1:3, x="topleft")
#' # See which points affect (5.5, 3 themost)
#' data.frame(X, cov=gp$kernel$k(X, c(5.5,3))) %>% arrange(-cov)
#' plot(k2b)
# LatentFactorKernel ----
LatentFactorKernel <- R6::R6Class(
  classname = "GauPro_kernel_LatentFactorKernel",
  inherit = GauPro_kernel,
  public = list(
    p = NULL, # vector
    p_est = NULL,
    p_lower = NULL,
    p_upper = NULL,
    p_length = NULL,
    s2 = NULL, # variance coefficient to scale correlation matrix to covariance
    s2_est = NULL,
    logs2 = NULL,
    logs2_lower = NULL,
    logs2_upper = NULL,
    nlevels = NULL,
    latentdim = NULL,
    xindex = NULL,
    pf_to_p_log = NULL,
    p_to_pf_inds = NULL,
    offdiagequal = NULL,
    #' @description Initialize kernel object
    #' @param p Vector of latent variables
    #' @param s2 Initial variance
    #' @param D Number of input dimensions of data
    #' @param p_lower Lower bound for p
    #' @param p_upper Upper bound for p
    #' @param p_est Should p be estimated?
    #' @param s2_lower Lower bound for s2
    #' @param s2_upper Upper bound for s2
    #' @param s2_est Should s2 be estimated?
    #' @param xindex Index of X to use the kernel on
    #' @param nlevels Number of levels for the factor
    #' @param latentdim Dimension of embedding space
    #' @param useC Should C code used? Much faster.
    #' @param offdiagequal What should offdiagonal values be set to when the
    #' indices are the same? Use to avoid decomposition errors, similar to
    #' adding a nugget.
    initialize = function(s2=1, D, nlevels, xindex,
                          latentdim,
                          p_lower=0, p_upper=1, p_est=TRUE,
                          s2_lower=1e-8, s2_upper=1e8, s2_est=TRUE,
                          useC=TRUE, offdiagequal=1-1e-6
    ) {
      # Must give in D
      if (missing(D)) {stop("Must give Index kernel D")}

      # latentdim defaults to 1 for D<=3, 2 for D>=4
      if (missing(latentdim)) {
        latentdim <- ifelse(D>3.5, 2, 1)
      }

      stopifnot(length(D) == 1, length(nlevels) == 1,
                length(xindex) == 1, length(latentdim) == 1,
                D>=1L, nlevels>=2L, xindex>=1L, latentdim>=1)
      # Following avoids redundancies
      stopifnot(latentdim < nlevels)

      self$D <- D
      self$nlevels <- nlevels
      self$xindex <- xindex
      self$latentdim <- latentdim

      p_to_pf <- c()
      for (i in 1:nlevels) {
        n0 <- max(0, latentdim+1-i)
        nnon0 <- latentdim - n0
        p_to_pf <- c(p_to_pf, rep(T, nnon0), rep(F, n0))
      }
      pf_to_p <- which(p_to_pf)
      # p_to_pf
      # pf_to_p
      # (1:(nlev*nld))[pf_to_p]
      # Names are backwards, or just unclear
      self$p_to_pf_inds <- pf_to_p
      self$pf_to_p_log  <- p_to_pf
      self$useC <- useC
      # Latent vars for first dim will be pinned to 0
      # p <- rnorm(latentdim*(nlevels - 1))
      # Now pinning more to 0, more complex conversions
      p <- rnorm(length(self$p_to_pf_inds))
      self$p <- p
      self$p_length <- length(p)
      # Ensure separation between levels to avoid instability
      self$p_lower <- rep(-25, self$p_length)
      # Don't give upper 1 since it will give optimization error
      self$p_upper <- rep(25, self$p_length)

      self$p_est <- p_est
      self$s2 <- s2
      self$logs2 <- log(s2, 10)
      self$logs2_lower <- log(s2_lower, 10)
      self$logs2_upper <- log(s2_upper, 10)
      self$s2_est <- s2_est

      self$offdiagequal <- offdiagequal
    },
    #' @description Calculate covariance between two points
    #' @param x vector.
    #' @param y vector, optional. If excluded, find correlation
    #' of x with itself.
    #' @param p Correlation parameters.
    #' @param s2 Variance parameter.
    #' @param params parameters to use instead of beta and s2.
    k = function(x, y=NULL, p=self$p, s2=self$s2, params=NULL) {
      if (!is.null(params)) {
        lenparams <- length(params)

        if (self$p_est) {
          p <- params[1:self$p_length]
        } else {
          p <- self$p
        }
        if (self$s2_est) {
          logs2 <- params[lenparams]
        } else {
          logs2 <- self$logs2
        }

        s2 <- 10^logs2
      } else {
        if (is.null(p)) {p <- self$p}
        if (is.null(s2)) {s2 <- self$s2}
      }

      pf <- self$p_to_pf(p)
      if (is.null(y)) {
        if (is.matrix(x)) {
          if (self$useC) {
            val <- s2 * corr_latentfactor_matrix_symC(x, pf, self$xindex,
                                                      self$latentdim,
                                                      self$offdiagequal)
          } else {
            val <- outer(1:nrow(x), 1:nrow(x),
                         Vectorize(function(i,j){
                           self$kone(x[i,],x[j,],pf=pf, s2=s2, isdiag=i==j)
                         }))
          }
          return(val)
        } else {
          return(s2 * 1)
        }
      }
      if (is.matrix(x) & is.matrix(y)) {
        # C took 0.000 sec, R took 1.793 sec
        if (self$useC) { # Way faster
          s2 * corr_latentfactor_matrixmatrixC(
            x=x, y=y, theta=pf, xindex=self$xindex,
            latentdim = self$latentdim, offdiagequal=self$offdiagequal)
        } else {
          outer(1:nrow(x), 1:nrow(y),
                Vectorize(function(i,j){self$kone(x[i,],y[j,],
                                                  pf=pf, s2=s2, isdiag=FALSE)}))
        }
      } else if (is.matrix(x) & !is.matrix(y)) {
        apply(x, 1, function(xx) {self$kone(xx, y, pf=pf, s2=s2)})
      } else if (is.matrix(y)) {
        apply(y, 1, function(yy) {self$kone(yy, x, pf=pf, s2=s2)})
      } else {
        self$kone(x, y, pf=pf, s2=s2)
      }
    },
    #' @description Find covariance of two points
    #' @param x vector
    #' @param y vector
    #' @param pf correlation parameters on regular scale, includes zeroes
    #' for first level.
    #' @param s2 Variance parameter
    #' @param isdiag Is this on the diagonal of the covariance?
    #' @param offdiagequal What should offdiagonal values be set to when the
    #' indices are the same? Use to avoid decomposition errors, similar to
    #' adding a nugget.
    #' @references https://stackoverflow.com/questions/27086195/linear-index-upper-triangular-matrix
    kone = function(x, y, pf, s2, isdiag=TRUE, offdiagequal=self$offdiagequal) {
      x <- x[self$xindex]
      y <- y[self$xindex]
      stopifnot(x>=1, y>=1, x<=self$nlevels, y<=self$nlevels,
                length(pf) == self$nlevels*self$latentdim,
                abs(x-as.integer(x)) < 1e-8, abs(y-as.integer(y)) < 1e-8)
      if (x==y) {
        # out <- s2 * 1
        # Trying to avoid singular values
        if (isdiag) {
          out <- s2 * 1
        } else {
          out <- s2 * offdiagequal
        }
      } else {
        # i <- x-1
        # j <- y-1
        # i <- min(x,y) #min(x-1, y-1)
        # j <- max(x,y) - 1 #max(x-1, y-1)
        # n <- self$nlevels
        # p_dist <- sum(p[i:j])
        latentx <- pf[(x-1)*self$latentdim+1:self$latentdim]
        latenty <- pf[(y-1)*self$latentdim+1:self$latentdim]
        p_dist2 <- sum((latentx - latenty)^2)
        out <- s2 * exp(-p_dist2)
      }
      if (any(is.nan(out))) {stop("Error #1341982")}
      out
    },
    #' @description Derivative of covariance with respect to parameters
    #' @param params Kernel parameters
    #' @param X matrix of points in rows
    #' @param C_nonug Covariance without nugget added to diagonal
    #' @param C Covariance with nugget
    #' @param nug Value of nugget
    dC_dparams = function(params=NULL, X, C_nonug, C, nug) {
      n <- nrow(X)

      stopifnot(X[, self$xindex] >= 1, X[, self$xindex] <= self$nlevels)

      lenparams <- length(params)

      if (lenparams > 0) {
        if (self$p_est) {
          p <- params[1:self$p_length]
        } else {
          p <- self$p
        }
        if (self$s2_est) {
          logs2 <- params[lenparams]
        } else {
          logs2 <- self$logs2
        }
      } else {
        p <- self$p
        logs2 <- self$logs2
      }

      log10 <- log(10)
      s2 <- 10 ^ logs2

      if (missing(C_nonug)) { # Assume C missing too, must have nug
        C_nonug <- self$k(x=X, params=params)
        C <- C_nonug + diag(nug*s2, nrow(C_nonug))
      }

      lenparams_D <- self$p_length*self$p_est + self$s2_est

      pf <- self$p_to_pf(p)

      if (self$useC) {
        dC_dparams <- kernel_latentFactor_dC(X, pf, C_nonug, self$s2_est,
                                             self$p_est, lenparams_D, s2*nug,
                                             self$latentdim, self$xindex-1,
                                             self$nlevels, s2)
      } else {
        dC_dparams <- array(dim=c(lenparams_D, n, n), data=0)
        if (self$s2_est) {
          dC_dparams[lenparams_D,,] <- C * log10
        }

        # Repeatedly calling self$ attributes is slow,
        # it's faster to just store as new variable
        latentdim <- self$latentdim
        xindex <- self$xindex

        if (self$p_est) {
          stopifnot(self$nlevels>=2L)
          for (k in 2:self$nlevels) { # k is index of level
            # kinds <- (k-1)*latentdim+1:latentdim - latentdim
            kinds <- (cumsum(self$pf_to_p_log) * self$pf_to_p_log)[
              (k-1)*latentdim+1:latentdim]
            kinds <- kinds[kinds != 0]
            stopifnot(length(kinds)>0, !anyDuplicated(kinds))
            # kactiveinds <- kinds[self$pf_to_p_log[kinds]]
            # kactiveinds <- ((kinds - 1) %% self$nlevels) + 1
            kactiveinds <- 1:length(kinds)
            stopifnot(length(kinds) == length(kactiveinds))
            for (i in seq(1, n-1, 1)) { # Index of X
              xlev <- X[i, xindex]
              latentx <- pf[(xlev-1)*latentdim+1:latentdim]
              for (j in seq(i+1, n, 1)) { # Index of Y
                ylev <- X[j, xindex]
                if (xlev > 1.5 && xlev == k && ylev != k) {
                  latenty <- pf[(ylev-1)*latentdim+1:latentdim]
                  p_dist2 <- sum((latentx - latenty)^2)
                  out <- s2 * exp(-p_dist2)
                  # kinds <- (xlev-1)*latentdim+1:latentdim - latentdim
                  # dC_dparams[kinds,i,j] <- -2 * out * (latentx - latenty)
                  # dC_dparams[kinds,j,i] <- dC_dparams[kinds,i,j]
                  dC_dparams[kinds,i,j] <- -2 * out * (latentx[kactiveinds] -
                                                         latenty[kactiveinds])
                  dC_dparams[kinds,j,i] <- dC_dparams[kinds,i,j]
                } else if (ylev > 1.5 && xlev != k && ylev == k) {
                  # latentx <- pf[(xlev-1)*latentdim+1:latentdim]
                  latenty <- pf[(ylev-1)*latentdim+1:latentdim]
                  p_dist2 <- sum((latentx - latenty)^2)
                  out <- s2 * exp(-p_dist2)
                  # kinds <- (ylev-1)*latentdim+1:latentdim - latentdim
                  # dC_dparams[kinds,i,j] <- 2 * out * (latentx - latenty)
                  # dC_dparams[kinds,j,i] <- dC_dparams[kinds,i,j]
                  dC_dparams[kinds,i,j] <- 2 * out * (latentx[kactiveinds] -
                                                        latenty[kactiveinds])
                  dC_dparams[kinds,j,i] <- dC_dparams[kinds,i,j]
                } else {
                  # Derivative is when when level isn't used in either
                  #  or when used in both.
                }
              }
            }
            for (i in seq(1, n, 1)) { # Get diagonal set to zero
              dC_dparams[k-1,i,i] <- 0
            }
          }
        }
      }
      return(dC_dparams)
    },
    #' @description Calculate covariance matrix and its derivative
    #'  with respect to parameters
    #' @param params Kernel parameters
    #' @param X matrix of points in rows
    #' @param nug Value of nugget
    C_dC_dparams = function(params=NULL, X, nug) {
      s2 <- self$s2_from_params(params)
      C_nonug <- self$k(x=X, params=params)
      C <- C_nonug + diag(s2*nug, nrow(X))
      dC_dparams <- self$dC_dparams(params=params, X=X, C_nonug=C_nonug, C=C, nug=nug)
      list(C=C, dC_dparams=dC_dparams)
    },
    #' @description Derivative of covariance with respect to X
    #' @param XX matrix of points
    #' @param X matrix of points to take derivative with respect to
    #' @param ... Additional args, not used
    dC_dx = function(XX, X, ...) {
      if (!is.matrix(XX)) {stop()}
      d <- ncol(XX)
      if (ncol(X) != d) {stop()}
      n <- nrow(X)
      nn <- nrow(XX)
      dC_dx <- array(0, dim=c(nn, d, n))
      dC_dx[, self$xindex, ] <- NA
      dC_dx
    },
    #' @description Starting point for parameters for optimization
    #' @param jitter Should there be a jitter?
    #' @param y Output
    #' @param p_est Is p being estimated?
    #' @param s2_est Is s2 being estimated?
    param_optim_start = function(jitter=F, y, p_est=self$p_est,
                                 s2_est=self$s2_est) {
      if (p_est) {
        vec <- pmin(pmax(self$p + jitter*rnorm(length(self$p), 0, 1),
                         self$p_lower), self$p_upper)
      } else {
        vec <- c()
      }
      if (s2_est) {
        vec <- c(vec, self$logs2 + jitter*rnorm(1))
      }
      vec
    },
    #' @description Starting point for parameters for optimization
    #' @param jitter Should there be a jitter?
    #' @param y Output
    #' @param p_est Is p being estimated?
    #' @param s2_est Is s2 being estimated?
    param_optim_start0 = function(jitter=F, y, p_est=self$p_est,
                                  s2_est=self$s2_est) {
      if (p_est) {
        vec <- pmin(pmax(jitter*rnorm(length(self$p), 0, 1),
                         self$p_lower), self$p_upper)
      } else {
        vec <- c()
      }
      if (s2_est) {
        vec <- c(vec, self$logs2 + jitter*rnorm(1))
      }
      vec
    },
    #' @description Lower bounds of parameters for optimization
    #' @param p_est Is p being estimated?
    #' @param s2_est Is s2 being estimated?
    param_optim_lower = function(p_est=self$p_est,
                                 s2_est=self$s2_est) {
      if (p_est) {vec <- c(self$p_lower)} else {vec <- c()}
      if (s2_est) {vec <- c(vec, self$logs2_lower)} else {}
      vec
    },
    #' @description Upper bounds of parameters for optimization
    #' @param p_est Is p being estimated?
    #' @param s2_est Is s2 being estimated?
    param_optim_upper = function(p_est=self$p_est,
                                 s2_est=self$s2_est) {
      if (p_est) {vec <- c(self$p_upper)} else {vec <- c()}
      if (s2_est) {vec <- c(vec, self$logs2_upper)} else {}
      vec
    },
    #' @description Set parameters from optimization output
    #' @param optim_out Output from optimization
    #' @param p_est Is p being estimated?
    #' @param s2_est Is s2 being estimated?
    set_params_from_optim = function(optim_out, p_est=self$p_est,
                                     s2_est=self$s2_est) {
      loo <- length(optim_out)
      if (p_est) {
        self$p <- optim_out[1:(self$p_length)]
      }
      if (s2_est) {
        self$logs2 <- optim_out[loo]
        self$s2 <- 10 ^ self$logs2
      }
    },
    #' @description Convert p (short parameter vector) to pf (long parameter
    #' vector with zeros).
    #' @param p Parameter vector
    p_to_pf = function(p) {
      pf <- rep(0, length(self$pf_to_p_log))
      pf[self$pf_to_p_log] <- p
      pf
    },
    #' @description Get s2 from params vector
    #' @param params parameter vector
    #' @param s2_est Is s2 being estimated?
    s2_from_params = function(params, s2_est=self$s2_est) {
      if (s2_est && !is.null(params)) { # Is last if in params
        10 ^ params[length(params)]
      } else { # Else it is just using set value, not being estimated
        self$s2
      }
    },
    #' @description Plot the points in the latent space
    plotLatent = function() {
      pf <- self$p_to_pf(self$p)
      pmat <- matrix(pf, ncol=self$latentdim, byrow=TRUE)
      pdf <- as.data.frame(pmat)
      pdf$name <- paste0("x=",1:nrow(pdf))
      if (self$latentdim == 1) {
        ggplot2::ggplot(pdf, ggplot2::aes(V1, 0, label=name)) +
          ggplot2::geom_point() +
          ggplot2::scale_y_continuous(breaks=NULL) +
          ggrepel::geom_label_repel() +
          ggplot2::ylab(NULL)
      } else if (self$latentdim == 2) {
        ggplot2::ggplot(pdf, ggplot2::aes(V1, V2, label=name)) +
          ggplot2::geom_point() +
          ggrepel::geom_label_repel()
      } else {
        stop("Can't plotLatent for latentdim > 2")
      }
    },
    #' @description Print this object
    print = function() {
      cat('GauPro kernel: Latent factor\n')
      cat('\tD  =', self$D, '\n')
      cat('\ts2 =', self$s2, '\n')
      cat('\ton x-index', self$xindex, 'with', self$nlevels, 'levels\n')
      cat('\t in', self$latentdim, 'latent dimensions\n')
    }
  )
)

