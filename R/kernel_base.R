# Kernels should implement:
# k kernel function for two vectors
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



#' Kernel R6 class
#'
#' @docType class
#' @importFrom R6 R6Class
# @export
#' @useDynLib GauPro, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats optim
# @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @field D Number of input dimensions of data
#' @field useC Should C code be used when possible? Can be much faster.
#' @examples
#' #k <- GauPro_kernel$new()
GauPro_kernel <- R6::R6Class(
  classname = "GauPro_kernel",
  public = list(
    D = NULL,
    useC = TRUE
    # k_diag = function(x) {
    #   if (is.matrix(x)) {rep(self$s2, nrow(x))}
    #   else if (length(x) == self$D) {self$s2}
    #   else {stop("Error in k_diag #4928")}
    # }
    ,
    # #' @description Plot kernel decay.
    # plot = function() {
    #   stopifnot(!is.null(self$D), self$D >= 1)
    #   n <- 51
    #   xseq <- seq(0,1,l=n)
    #   x0 <- rep(0, self$D)
    #   df <- NULL
    #   for (i in 1:self$D) {
    #     X <- matrix(0, ncol=self$D, nrow=n)
    #     X[, i] <- xseq
    #     # xi <- rep(0, self$D)
    #     # xi[i] <-
    #     k <- self$k(x0, X)
    #     df <- rbind(df,
    #                 data.frame(i=i, x2=xseq, k=k)
    #     )
    #   }
    #   ggplot2::ggplot(df, ggplot2::aes(x2, k)) + ggplot2::geom_line() +
    #     ggplot2::facet_wrap(.~i)
    # },
    #' @description Plot kernel decay.
    #' @param X Matrix of points the kernel is used with. Some will be used
    #' to demonstrate how the covariance changes.
    plot = function(X=NULL) {
      stopifnot(!is.null(self$D), self$D >= 1)
      if (!is.null(X)) {
        stopifnot(is.matrix(X), ncol(X)==self$D)
      }
      n <- 101
      x0 <- rep(1, self$D)
      # df <- NULL
      plots <- list()
      factorinfo <- GauPro:::find_kernel_factor_dims(self)
      factordims <- if (length(factorinfo)==0) {numeric()} else {
        factorinfo[seq(1, length(factorinfo), 2)]}
      # Loop over each dimension
      for (i in 1:self$D) {
        df <- NULL
        if (!(i %in% factordims)) {
          if (is.null(X)) {
            Xi <- seq(0, 1, l=10)
          } else {
            Xi <- X[, i]
          }
          minXi <- min(Xi)
          maxXi <- max(Xi)
          Xiseq <- seq(minXi, maxXi,l=n)
          XX <- matrix(rep(x0, n), byrow = T, ncol=length(x0))
          # u is values to use as centers for current dimension
          u <- seq(minXi, maxXi, l=3)
          for (j in seq_along(u)) {
            x0j <- x0
            x0j[i] <- u[j]
            XX[, i] <- Xiseq
            k <- self$k(x0j, XX)
            df <- rbind(df,
                        data.frame(i=i, x1i=j,
                                   x1=u[j], x2=Xiseq, k=k)
            )
          }
          ploti <- ggplot2::ggplot(df,
                                   ggplot2::aes(x2, k,
                                                group=x1, color=factor(x1i))) +
            ggplot2::geom_line() +
            # ggplot2::facet_wrap(.~i, scales='free_x') +
            ggplot2::guides(color='none') +
            ggplot2::xlab(paste0("X", i)) +
            ggplot2::coord_cartesian(ylim=(c(0, max(df$k))))
        } else { # Factor dim
          # Copied from kernel_LatentFactor plot
          # nlevels <- length(unique(X[, i]))
          nlevels <- which(factorinfo[seq(1, length(factorinfo), 2)] == i)[1] *2
          xindex <- i
          x1 <- 1:nlevels
          X1 <- X2 <- matrix(data=1, ncol=self$D, nrow=nlevels)
          X1[, xindex] <- x1
          X2[, xindex] <- x1
          k <- self$k(X1, X2)

          df <- NULL
          for (il in 1:nlevels) {
            for (j in 1:nlevels) {
              df <- rbind(df,
                          data.frame(x1=il, x2=j, k=k[il,j]))
            }
          }
          ploti <- ggplot2::ggplot(data=df, ggplot2::aes(x1, x2, fill=k)) +
            ggplot2::geom_tile() +
            ggplot2::scale_fill_gradient(low='white', high='red', limits=c(0, NA)) +
            ggplot2::scale_y_reverse() +
            ggplot2::ylab(NULL) +
            ggplot2::xlab(paste0("X", i))

        }
        # Add plot to list
        plots[[i]] <- ploti
        rm(df)
      }
      # Arrange
      # do.call(gridExtra::grid.arrange, c(plots, ncol=floor(sqrt(length(plots)))))
      gridExtra::grid.arrange(grobs=plots,
                              ncol=floor(sqrt(length(plots))))
    },
    # plot = function(X=NULL) {
    #   stopifnot(!is.null(self$D), self$D >= 1)
    #   if (!is.null(X)) {
    #     stopifnot(is.matrix(X), ncol(X)==self$D)
    #   }
    #   n <- 101
    #   x0 <- rep(0, self$D)
    #   df <- NULL
    #   # Loop over each dimension
    #   for (i in 1:self$D) {
    #     if (is.null(X)) {
    #       Xi <- seq(0, 1, l=10)
    #     } else {
    #       Xi <- X[, i]
    #     }
    #     minXi <- min(Xi)
    #     maxXi <- max(Xi)
    #     Xiseq <- seq(minXi, maxXi,l=n)
    #     XX <- matrix(rep(x0, n), byrow = T, ncol=length(x0))
    #     # u is values to use as centers for current dimension
    #     u <- seq(minXi, maxXi, l=3)
    #     for (j in seq_along(u)) {
    #       x0j <- x0
    #       x0j[i] <- u[j]
    #       XX[, i] <- Xiseq
    #       k <- self$k(x0j, XX)
    #       df <- rbind(df,
    #                   data.frame(i=i, x1i=j,
    #                              x1=u[j], x2=Xiseq, k=k)
    #       )
    #     }
    #   }
    #   ggplot2::ggplot(df, ggplot2::aes(x2, k, group=x1, color=factor(x1i))) +
    #     ggplot2::geom_line() +
    #     ggplot2::facet_wrap(.~i, scales='free_x') +
    #     ggplot2::guides(color='none')
    # },
    #' @description Print this object
    print = function() {
      cat('GauPro kernel: (type unknown)\n')
      cat('\tD =', self$D, '\n')
    }
  ),
  private = list(

  )
)
