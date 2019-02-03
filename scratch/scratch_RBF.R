# RBFcreate <- function(X, Y) {
#   self <- list()
#   self$X
#   self$Y
#   self
# }

RBFfit <- function(X, Y) {#browser()
  self <- list()
  centers <- X
  n <- nrow(X)
  self$BasisFunc <- function(x,y) {
    sqrt(abs((x - y)^2))
  }
  Phi <- outer(1:n, 1:n, Vectorize(
    function(i,j) {
      # sqrt(abs((X[i,] - X[j,])^2))
      self$BasisFunc(X[i,], X[j,])
    }
  ))
  Phiinv_Y <- solve(Phi, Y)

  self$X <- X
  self$Y <- Y
  self$Phiinv_Y <- Phiinv_Y
  self$Phi <- Phi
  self
}

RBFpred <- function(self, xp) {#browser()
  n <- nrow(self$X)
  np <- nrow(xp)
  phi <- outer(1:n, 1:np, function(i, ip) {
    # sqrt(abs((self$X[i,] - xp[ip,])^2))
    self$BasisFunc(self$X[i,], xp[ip,])
  })
  Ypredmean <- t(phi) %*% self$Phiinv_Y

  # Get var
  Ypredvar <- self$BasisFunc(self$X[1,], self$X[1,]) - t(phi) %*% solve(self$Phi, phi)

  # Return mean and var
  list(mean=Ypredmean, var=diag(-Ypredvar), cov=Ypredvar)
}


d <- 1
n <- 10
x <- matrix(runif(d*n), ncol=d)
f <- function(x) {cos(x^1.9)}
y <- as.matrix(apply(x, 1, f))
rb1 <- RBFfit(X=x, Y = y)
np <- 30
xp <- matrix(runif(d*np), ncol=d)
RBFpred(rb1, xp)

xl <- seq(0,1,l=100)
pr <- RBFpred(rb1, as.matrix(xl))
plot(xl, pr$me+2*sqrt(pr$v), type='l', col=3, ylim=c(min(pr$me-2*sqrt(pr$v)), max(pr$me+2*sqrt(pr$v))))
points(xl, pr$me-2*sqrt(pr$v), type='l', col=3)
points(xl, pr$me, type='l')
points(x, y)









AugRBFfit <- function(X, Y) {#browser()
  self <- list()
  centers <- X
  n <- nrow(X)
  self$BasisFunc <- function(x,y) {
    sqrt(abs((x - y)^2))
  }
  Phi <- outer(1:n, 1:n, Vectorize(
    function(i,j) {
      # sqrt(abs((X[i,] - X[j,])^2))
      self$BasisFunc(X[i,], X[j,])
    }
  ))
  k <- 2
  P <- matrix(0, n, k*d+1)
  P[,1] <- 1
  # browser()
  for (ki in 1:k) {
    P[,1+ki*d - d + 1:d] <- X ^ ki
  }
  AugPhi <- rbind(cbind(Phi, P), cbind(t(P), matrix(0,k+1,k+1)))
  AugPhiinv_Y <- solve(AugPhi, rbind(Y, matrix(0,k*d+1,1)))
  theta <- AugPhiinv_Y[1:n]
  mu <- AugPhiinv_Y[(n+1):(n+k+1)]

  self$X <- X
  self$Y <- Y
  self$k <- k
  self$AugPhiinv_Y <- AugPhiinv_Y
  self$Phi <- Phi
  self$theta <- theta
  self$mu <- mu
  self
}

AugRBFpred <- function(self, xp) {#browser()
  print("using augmented")
  n <- nrow(self$X)
  np <- nrow(xp)
  phi <- outer(1:n, 1:np, function(i, ip) {
    # sqrt(abs((self$X[i,] - xp[ip,])^2))
    self$BasisFunc(self$X[i,], xp[ip,])
  })
  # browser()
  Ypredmean <- t(phi) %*% self$theta + rowSums(matrix(self$mu, nrow(xp), self$k+1, byrow=T) * do.call(cbind, lapply(0:self$k, function(ki) xp^ki)))

  # Get var
  Ypredvar <- self$BasisFunc(self$X[1,], self$X[1,]) - t(phi) %*% solve(self$Phi, phi)

  # Return mean and var
  list(mean=Ypredmean, var=diag(-Ypredvar), cov=Ypredvar)
}



d <- 1
n <- 10
x <- matrix(runif(d*n), ncol=d)
f <- function(x) {cos(x^1.9)}
y <- as.matrix(apply(x, 1, f))
arb1 <- AugRBFfit(X=x, Y = y)
np <- 30
xp <- matrix(runif(d*np), ncol=d)
AugRBFpred(arb1, xp)

xl <- seq(0,1,l=100)
apr <- AugRBFpred(arb1, as.matrix(xl))
plot(xl, apr$me+2*sqrt(apr$v), type='l', col=3, ylim=c(min(apr$me-2*sqrt(apr$v)), max(apr$me+2*sqrt(apr$v))))
points(xl, apr$me-2*sqrt(apr$v), type='l', col=3)
points(xl, apr$me, type='l')
points(x, y)
