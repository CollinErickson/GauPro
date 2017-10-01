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
  if (x$D == 1) {
    x$cool1Dplot(...)
  } else if (x$D == 2) {
    x$plot2D(...)
  } else {
    stop("No plot method for higher than 2 dimension")
  }
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
  kernel_product$new(k1=k1, k2=k2)
}
