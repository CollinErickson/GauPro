#' Predict for class GauPro
#'
#' @param gp Object of class GauPro
#' @param XX new points to predict
#' @param se.fit Should standard error be returned (and variance)?
#' @param covmat Should the covariance matrix be returned?
#' @param split_speed Should the calculation be split up to speed it up?
#'
#' @return
#' @export
#'
#' @examples
#' n <- 12
#' x <- matrix(seq(0,1,length.out = n), ncol=1)
#' y <- sin(2*pi*x) + rnorm(n,0,1e-1)
#' gp <- GauPro(X=x, Z=y, parallel=FALSE)
#' predict(gp, .448)
predict.GauPro <- function(gp, XX, se.fit=F, covmat=F, split_speed=T) {
  gp$predict(XX=XX, se.fit=se.fit, covmat=covmat, split_speed=split_speed)
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
#' k1 <- Gaussian$new(theta=1)
#' k2 <- Gaussian$new(theta=2)
#' k <- k1 + k2
'+.GauPro_kernel' <- function(k1, k2) {
  kernel_sum$new(k1=k1, k2=k2)
}
