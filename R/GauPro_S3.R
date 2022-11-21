# S3 methods for GauPro_kernel_model, which has class GauPro
# plot, print, and format are automatically dispatched, all others must be added

#' Predict for class GauPro
#'
#' @param object Object of class GauPro
#' @param XX new points to predict
#' @param se.fit Should standard error be returned (and variance)?
#' @param covmat Should the covariance matrix be returned?
#' @param split_speed Should the calculation be split up to speed it up?
#' @param ... Additional parameters
#'
#' @return Prediction from object at XX
#' @export
#'
#' @examples
#' n <- 12
#' x <- matrix(seq(0,1,length.out = n), ncol=1)
#' y <- sin(2*pi*x) + rnorm(n,0,1e-1)
#' gp <- GauPro(X=x, Z=y, parallel=FALSE)
#' predict(gp, .448)
predict.GauPro <- function(object, XX, se.fit=F, covmat=F, split_speed=T, ...) {
  object$predict(XX=XX, se.fit=se.fit, covmat=covmat, split_speed=split_speed)
}

if (F) {
  # Plot is automatically dispatched, same with print and format
  #' Plot for class GauPro
  #'
  #' @param x Object of class GauPro
  #' @param ... Additional parameters
  #'
  #' @return Nothing
  #' @export
  #'
  #' @examples
  #' n <- 12
  #' x <- matrix(seq(0,1,length.out = n), ncol=1)
  #' y <- sin(2*pi*x) + rnorm(n,0,1e-1)
  #' gp <- GauPro(X=x, Z=y, parallel=FALSE)
  #' if (requireNamespace("MASS", quietly = TRUE)) {
  #'   plot(gp)
  #' }
  #'
  plot.GauPro <- function(x,  ...) {
    x$plot(...)
    # if (x$D == 1) {
    #   x$cool1Dplot(...)
    # } else if (x$D == 2) {
    #   x$plot2D(...)
    # } else {
    #   # stop("No plot method for higher than 2 dimension")
    #   x$plotmarginal()
    # }
  }
}

#' Summary for GauPro object
#'
#' @param object GauPro R6 object
#' @param ... Additional arguments passed to summary
#'
#' @return Summary
#' @export
summary.GauPro <- function(object, ...) {
  object$summary(...)
}

#' Print summary.GauPro
#'
#' @param x summary.GauPro object
#' @param ... Additional args
#' @importFrom stats binom.test
#'
#' @return
#' @export
print.summary.GauPro <- function(x, ...) {
  # Formula
  cat("Formula:\n")
  cat("\t", x$formula, "\n\n")

  # Residuals
  cat("Residuals:\n")
  print(summary(x$residualsLOO))

  # Importance
  cat("\nFeature importance:\n")
  print(x$importance)

  # R-squared, Adj R-squared
  cat("\nPseudo leave-one-out R-squared:\n")
  cat("\t", x$r.squaredLOO, "\n")

  # Coverage
  # cat("\nLeave-one-out 95% coverage:\n")
  # cat("\t", x$coverageLOO, "\t(on", x$N, "samples)", "\n")
  pval68 <- binom.test(x$coverage68LOO*x$N, x$N, .68)$p.value
  pval95 <- binom.test(x$coverage95LOO*x$N, x$N, .95)$p.value
  cat("\nLeave-one-out coverage (on", x$N, "samples):\n")
  cat("\t68%:  ", x$coverage68LOO, "\t\tp-value:  ", pval68, "\n")
  cat("\t95%:  ", x$coverage95LOO, "\t\tp-value:  ", pval95, "\n")

  # Return invisible self
  invisible(x)
}


#' Kernel sum
#'
#' @param k1 First kernel
#' @param k2 Second kernel
#'
#' @return Kernel which is sum of two kernels
#' @export
#'
#' @examples
#' k1 <- Exponential$new(beta=1)
#' k2 <- Matern32$new(beta=0)
#' k <- k1 + k2
#' k$k(matrix(c(2,1), ncol=1))
'+.GauPro_kernel' <- function(k1, k2) {
  if (is.numeric(k1) && k1==0) {
    return(k2)
  }
  if (is.numeric(k2) && k2==0) {
    return(k1)
  }
  if (!("GauPro_kernel" %in% class(k1))) {
    stop("Can only add GauPro kernels with other kernels")
  }
  if (!("GauPro_kernel" %in% class(k2))) {
    stop("Can only add GauPro kernels with other kernels")
  }
  kernel_sum$new(k1=k1, k2=k2)
}


#' Kernel product
#'
#' @param k1 First kernel
#' @param k2 Second kernel
#'
#' @return Kernel which is product of two kernels
#' @export
#'
#' @examples
#' k1 <- Exponential$new(beta=1)
#' k2 <- Matern32$new(beta=0)
#' k <- k1 * k2
#' k$k(matrix(c(2,1), ncol=1))
'*.GauPro_kernel' <- function(k1, k2) {
  if (is.numeric(k1) && k1==1) {
    return(k2)
  }
  if (is.numeric(k2) && k2==1) {
    return(k1)
  }
  if (!("GauPro_kernel" %in% class(k1))) {
    stop("Can only multiply GauPro kernels with other kernels")
  }
  if (!("GauPro_kernel" %in% class(k2))) {
    stop("Can only multiply GauPro kernels with other kernels")
  }
  kernel_product$new(k1=k1, k2=k2)
}
