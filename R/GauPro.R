library(R6)
GauPro <- R6Class(classname = "GauPro",
  public = list(
    X = NA,
    Z = NA,
    N = NA,
    D = NA,
    corr = "Gauss",
    nug = 1e-6,
    theta = 1,
    mu_hat = NA,
    s2_hat = NA,
    K = NA,
    Kchol = NA,
    Kinv = NA,
    fit = function(X, Z) {#browser()
      self$X <- X
      self$Z <- Z
      self$N <- nrow(X)
      self$D <- ncol(X)
      #self$update_params()
      self$theta_update()
    },
    update_params = function () {
      self$K <- gauss_cor_mat(self$X, theta=self$theta) + diag(self$nug, self$N)
      self$Kchol <- chol(self$K)
      self$Kinv <- chol2inv(self$Kchol)
      self$mu_hat <- sum(self$Kinv %*% self$Z) / sum(self$Kinv)
      self$s2_hat <- c(t(self$Z - self$mu_hat) %*% self$Kinv %*% (self$Z - self$mu_hat) / self$N)
    },
    pred = function(XX, covmat=F) {#browser()
      if (!is.matrix(XX)) {
        if (self$D == 1) XX <- matrix(XX, ncol=1)
        else if (length(XX) == self$D) XX <- matrix(XX, nrow=1)
        else stop('Predict input should be matrix')
      }
      #covmat <- gauss_cor(c(x, xx))
      #kx <- gauss_cor_mat(self$X) + diag(self$nug, self$N)
      kxx <- gauss_cor_mat(XX, theta=self$theta)
      kx.xx <- gauss_cor_mat(XX, self$X, theta=self$theta)

      mn <- self$pred_mean(XX, kx.xx=kx.xx)#self$mu_hat + kx.xx %*% self$Kinv %*% (self$Z - self$mu_hat)
      #s2 <- self$s2_hat * (1 - kx.xx %*% self$Kinv %*% t(kx.xx) + (1-sum(self$Kinv %*% t(kx.xx)))^2 / (sum(self$Kinv)))
      s2 <- self$pred_var(XX, kxx=kxx, kx.xx=kx.xx) #self$s2_hat * diag(kxx - kx.xx %*% self$Kinv %*% t(kx.xx))
      se <- sqrt(s2)
      if (covmat) return(list(mean=mn, s2=s2, se=se, cov=self$pred_var(XX, kxx=kxx, kx.xx=kx.xx, covmat=T)))
      data.frame(mean=mn, s2=s2, se=se)
    },
    pred_mean = function(XX, kx.xx) {
      self$mu_hat + kx.xx %*% self$Kinv %*% (self$Z - self$mu_hat)
    },
    pred_var = function(XX, kxx, kx.xx, covmat=F) {
      if (covmat) return(self$s2_hat * (kxx - kx.xx %*% self$Kinv %*% t(kx.xx)))
      self$s2_hat * diag(kxx - kx.xx %*% self$Kinv %*% t(kx.xx))
    },
    deviance = function (theta) {
      K <- gauss_cor_mat(self$X, theta=theta) + diag(self$nug, self$N)
      Kchol <- chol(K)
      Kinv <- chol2inv(Kchol)
      mu_hat <- sum(Kinv %*% self$Z) / sum(Kinv)
      s2_hat <- c(t(self$Z - self$mu_hat) %*% Kinv %*% (self$Z - mu_hat) / self$N)
      detK <- prod(diag(Kchol))^2
      log(detK) + self$N * log((self$Z - mu_hat) %*% solve(K, self$Z - mu_hat))
    },
    deviance_log = function (beta) {
      #points(beta,0)
      theta <- 10^beta
      self$deviance(theta=theta)
    },
    deviance_search = function () {
      #optimize(self$deviance_log, interval = c(-10,10), maximum = F, tol=1e-3)
      nvars <- self$D
      suppressWarnings(
        rgenoud::genoud(self$deviance_log, nvars=nvars,#self$D,#1,
                      Domains = matrix(c(-10,10), nvars, 2, byrow=T), max = F,
                      solution.tolerance=1e-3,
                      pop.size=min(1000, 200*self$D),
                      max.generations = min(10, 3*self$D),
                      print.level = 0
                      )
      )
    },
    theta_update = function () {
      dev <- self$deviance_search()
      #self$theta <- 10 ^ self$deviance_search()$minimum
      self$theta <- 10 ^ self$deviance_search()$par
      self$update_params()
    },
    cool1Dplot = function () {#browser()
      if (self$D != 1) stop('Must be 1D')
      minx <- min(self$X)
      maxx <- max(self$X)
      x1 <- minx - .1 * (maxx - minx)
      x2 <- maxx + .1 * (maxx - minx)
      nn <- 201
      x <- seq(x1, x2, length.out = nn)
      px <- self$pred(x, covmat = T)
      n2 <- 20
      newy <- MASS::mvrnorm(n=n2, mu=px$mean, Sigma=px$cov)
      plot(x,px$me, type='l', lwd=4, ylim=c(min(newy),max(newy)))
      sapply(1:n2, function(i) points(x, newy[i,], type='l', col='gray'))
      points(self$X, self$Z, pch=19, col=1, cex=2)
    }
  ),
  private = list(

  )
)
