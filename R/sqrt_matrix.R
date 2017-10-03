#' Find the square root of a matrix
#'
#' Same thing as 'expm::sqrtm', but faster.
#'
#' @param mat Matrix to find square root matrix of
#' @param symmetric Is it symmetric? Passed to eigen.
#'
#' @return Square root of mat
#' @export
#'
#' @examples
#' mat <- matrix(c(1,.1,.1,1), 2, 2)
#' smat <- sqrt_matrix(mat=mat, symmetric=TRUE)
#' smat %*% smat
sqrt_matrix = function(mat, symmetric) {
  e <- eigen(mat, symmetric=symmetric)
  V <- e$vectors
  if (length(V) == 1) { # diag in 1D is scalar and doesn't work correctly
    B <- sqrt(e$values) * V %*% t(V)
  } else {
    B <- V %*% diag(sqrt(e$values)) %*% t(V)
  }
  B
}
