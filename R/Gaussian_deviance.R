
Gaussian_devianceR = function (theta, nug, X, Z) { browser()# joint deviance
  N <- nrow(X)
  K <- corr_gauss_matrix_sym_armaC(X, theta=theta) + diag(nug, N)
  Kchol <- try(chol(K), silent = T)
  if (inherits(Kchol, "try-error")) {return(Inf)}
  Kinv <- chol2inv(Kchol)
  mu_hat <- sum(Kinv %*% Z) / sum(Kinv)
  logdetK <- 2 * sum(log(diag(Kchol)))
  logdetK + N * log(t(Z - mu_hat) %*% (Kinv %*% (Z - mu_hat)))
}
Gaussian_deviance_partC = function (theta, nug, X, Z) {browser()
  # Not faster than devianceC or even deviance
  N <- nrow(X)
  K <- corr_gauss_matrix_sym_armaC(X, theta=theta) + diag(nug, N)
  Kchol <- try(chol(K), silent = T)
  if (inherits(Kchol, "try-error")) {return(Inf)}
  Kinv <- chol2inv(Kchol)
  logdetK <- 2 * sum(log(diag(Kchol)))
  logdetK + Gaussian_deviance_part(theta, nug, X, Z, Kinv)
}
# This is fully C, see Gaussian_deviance.cpp file
#Gaussian_devianceC = function (theta, nug) {
  # 50% faster than deviance_partC and deviance on small,  twice on big
#  N <- nrow(X)
#  K <- corr_gauss_matrix_sym_armaC(X, theta=theta) + diag(nug, N)
#  devianceC(theta, nug, X, Z, K)
#}


Gaussian_deviance_gradR = function (theta, nug, X, Z) {browser()
  # Works but DONT USE since deviance_gradC is 10x faster, and only does theta grad so no good
  N <- nrow(X)
  D <- ncol(X)
  K <- corr_gauss_matrix_sym_armaC(X, theta=theta) + diag(nug, N)
  Kchol <- try(chol(K))
  if (inherits(Kchol, "try-error")) {return(Inf)}
  Kinv <- chol2inv(Kchol)
  mu_hat <- sum(Kinv %*% Z) / sum(Kinv)
  y <- Z - mu_hat
  Kinv.y <- Kinv %*% y
  t2a <- -N / (t(y)%*%Kinv.y)
  dD <- rep(NA,D)
  for(i in 1:D) {
    dK <- K
    for(j in 1:N) {
      for(k in 1:N) {
        dK[j, k] <- -(X[j,i]-X[k,i])^2 * dK[j, k]
      }
    }
    t1 <- sum(sapply(1:N, function(ii) {sum(Kinv[ii,] * dK[,ii])}))
    t2 <- t2a * t(Kinv.y) %*% dK %*% Kinv.y
    dD[i] <- 2 * (t1+t2)
  }
  dD
}
Gaussian_deviance_grad_fixR = function (theta, nug, X, Z) {
  # Grad doesn't include d(mu)/d(theta), tried adding it here
  #  If eqns are right, it makes a 1e-16 effect, so not worth it
  N <- nrow(X)
  D <- ncol(X)
  K <- corr_gauss_matrix_sym_armaC(X, theta=theta) + diag(nug, N)
  Kchol <- try(chol(K))
  if (inherits(Kchol, "try-error")) {return(Inf)}
  Kinv <- chol2inv(Kchol)
  mu_hat <- sum(Kinv %*% Z) / sum(Kinv)
  y <- Z - mu_hat
  Kinv.y <- Kinv %*% y
  t2a <- -N / (t(y)%*%Kinv.y)
  dD <- rep(NA,D)
  for(i in 1:D) {
    dK <- K
    for(j in 1:N) {
      for(k in 1:N) {
        dK[j, k] <- -(X[j,i]-X[k,i])^2 * dK[j, k]
      }
    }
    # mu = A^-1 %*% B, use it to calculate grad more accurately
    A <- sum(Kinv)
    B <- sum(Kinv.y)
    dA <- sum(Kinv)^2 * sum(Kinv %*% dK %*% Kinv)
    dB <- sum(Kinv %*% dK %*% Kinv %*%y)
    dmu = (A * dB - B * dA) / A^2
    t1 <- sum(sapply(1:N, function(ii) {sum(Kinv[ii,] * dK[,ii])}))
    t2 <- t2a * (t(Kinv.y) %*% dK %*% Kinv.y + - 2 * sum(Kinv.y) * dmu) # after including mu
    print(c(t(Kinv.y) %*% dK %*% Kinv.y , - 2 * sum(Kinv.y) * dmu))
    #t2 <- t2a * t(Kinv.y) %*% dK %*% Kinv.y # before mu
    dD[i] <- 2 * (t1+t2)
  }
  dD
}
Gaussian_deviance_gradC = function (theta, nug, X, Z, overwhat) {browser()
  N <- nrow(X)
  K <- corr_gauss_matrix_sym_armaC(X, theta=theta) + diag(nug, N)
  Kchol <- try(chol(K), silent = T)
  if (inherits(Kchol, "try-error")) {return(Inf)}
  Kinv <- chol2inv(Kchol)
  mu_hat <- sum(Kinv %*% Z) / sum(Kinv)
  y <- Z - mu_hat
  # Call RcppArmadillo function to calculate grad
  if (overwhat == "theta") {
    return(deviance_grad_theta(X, K, Kinv, y))
  } else if (overwhat == "nug") {
    return(deviance_grad_nug(X, K, Kinv, y))
  } else if (overwhat == "joint") {
    return(deviance_grad_joint(X, K, Kinv, y))
  } else {stop("No overwhat given #231854")}
}






Gaussian_deviance_fngrR = function (theta, nug, X, Z, overwhat) {
  N <- nrow(X)
  K <- corr_gauss_matrix_sym_armaC(X, theta=theta) + diag(nug, N)
  Kchol <- try(chol(K), silent = T)
  if (inherits(Kchol, "try-error")) {return(list(fn=Inf, gr=NA))}
  Kinv <- chol2inv(Kchol)
  mu_hat <- sum(Kinv %*% Z) / sum(Kinv)
  logdetK <- 2 * sum(log(diag(Kchol)))
  fn <- logdetK + N * log(t(Z - mu_hat) %*% (Kinv %*% (Z - mu_hat)))

  y <- Z - mu_hat
  # Call RcppArmadillo function to calculate grad
  if (overwhat == "theta") {
    gr <- (deviance_grad_theta(X, K, Kinv, y))
  } else if (overwhat == "nug") {
    gr <- (deviance_grad_nug(X, K, Kinv, y))
  } else if (overwhat == "joint") {
    gr <- (deviance_grad_joint(X, K, Kinv, y))
  } else {stop("No overwhat given #8523093")}

  return(list(fn=fn, gr=gr))
}
Gaussian_deviance_fngrC = function (theta, nug, X, Z, overwhat) {
  N <- nrow(X)
  K <- corr_gauss_matrix_sym_armaC(X, theta=theta) + diag(nug, N)
  # Call RcppArmadillo function to calculate all
  if (overwhat == "theta") {
    fngr <- (deviance_fngr_theta(X, Z, K))
  } else if (overwhat == "nug") {
    fngr <- (deviance_fngr_nug(X, Z, K))
  } else if (overwhat == "joint") {
    fngr <- (deviance_fngr_joint(X, Z, K))
  } else {stop("No overwhat given #935277")}

  return(list(fn=fngr[1], gr=fngr[-1]))
}
