#' GauPro_selector
#'
#' @param type Type of Gaussian process, or the kind of correlation function.
#' @param ... Pass on
#'
#' @return A GauPro object
#' @export
#'
#' @examples
#' n <- 12
#' x <- matrix(seq(0,1,length.out = n), ncol=1)
#' #y <- sin(2*pi*x) + rnorm(n,0,1e-1)
#' y <- (2*x) %%1
#' gp <- GauPro(X=x, Z=y)
GauPro <- function(..., type="Gauss") {
  gp <- GauPro_Gauss$new(...)
  return(gp)
}
