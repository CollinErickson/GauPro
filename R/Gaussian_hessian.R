#' Calculate Hessian for a GP with Gaussian correlation
#'
#' @param XX The vector at which to calculate the Hessian
#' @param X The input points
#' @param Z The output values
#' @param Kinv The inverse of the correlation matrix
#' @param mu_hat Estimate of mu
#' @param theta Theta parameters for the correlation
#'
#' @return Matrix, the Hessian at XX
#' @export
#'
#' @examples
#' set.seed(0)
#' n <- 40
#' x <- matrix(runif(n*2), ncol=2)
#' f1 <- function(a) {sin(2*pi*a[1]) + sin(6*pi*a[2])}
#' y <- apply(x,1,f1) + rnorm(n,0,.01)
#' gp <- GauPro(x,y, verbose=2, parallel=FALSE);gp$theta
#' gp$hessian(c(.2,.75), useC=FALSE) # Should be -38.3, -5.96, -5.96, -389.4 as 2x2 matrix
Gaussian_hessianR <- function(XX, X, Z, Kinv, mu_hat, theta) {#browser()
  n <- nrow(X) # number of points already in design
  d <- length(XX) # input dimensions
  Kinv_Zmu <- Kinv %*% (Z - mu_hat) #solve(R, Z - mu_hat)
  #d2K <- matrix()
  d2ZZ <- matrix(NA, d, d)
  #d2r <- numeric(n)
  exp_sum <- numeric(n)
  for (j in 1:n) {
    exp_sum[j] <- exp(-sum(theta * (XX - X[j, ])^2))
  }
  for (i in 1:d) { # diagonal points
    d2K_dxidxi <- numeric(n)
    for (j in 1:n) {

      d2Kj_dxidxi <- (-2 * theta[i] + 4 * theta[i]^2 * (XX[i] - X[j, i])^2) * exp_sum[j]
      d2K_dxidxi[j] <- d2Kj_dxidxi
    }

    tval <- t(d2K_dxidxi) %*% Kinv_Zmu
    d2ZZ[i, i] <- tval
  }
  if (d > 1) {
    for (i in 1:(d-1)) { # off diagonal points
      for (k in (i+1):d) {

        d2K_dxidxk <- numeric(n)
        for (j in 1:n) {

          d2Kj_dxidxk <- 4 * theta[i] * theta[k] * (XX[i] - X[j, i]) * (XX[k] - X[j, k]) * exp_sum[j]
          d2K_dxidxk[j] <- d2Kj_dxidxk
        }

        tval <- t(d2K_dxidxk) %*% Kinv_Zmu
        d2ZZ[i, k] <- tval
        d2ZZ[k, i] <- tval
      }
    }
  }


  hess <- d2ZZ #t(d2K) %*% solve(R, Z - mu_hat)
}


#' Calculate Hessian for a GP with Gaussian correlation
#'
#' @param XX The vector at which to calculate the Hessian
#' @param X The input points
#' @param Z The output values
#' @param Kinv The inverse of the correlation matrix
#' @param mu_hat Estimate of mu
#' @param theta Theta parameters for the correlation
#'
#' @return Matrix, the Hessian at XX
#' @export
#'
#' @examples
#' set.seed(0)
#' n <- 40
#' x <- matrix(runif(n*2), ncol=2)
#' f1 <- function(a) {sin(2*pi*a[1]) + sin(6*pi*a[2])}
#' y <- apply(x,1,f1) + rnorm(n,0,.01)
#' gp <- GauPro(x,y, verbose=2, parallel=FALSE);gp$theta
#' gp$hessian(c(.2,.75), useC=TRUE) # Should be -38.3, -5.96, -5.96, -389.4 as 2x2 matrix
Gaussian_hessianC <- function(XX, X, Z, Kinv, mu_hat, theta) {#browser()
  print("Using C version")
  #browser()
  Gaussian_hessianCC(XX, X, Z, Kinv, mu_hat, theta)
}
