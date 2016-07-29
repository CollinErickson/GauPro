library(R6)
GauPro <- R6Class(classname = "GauPro",
  public = list(
    X = NA,
    Z = NA,
    N = NA,
    D = NA,
    corr = "Gauss",
    nug = 1e-6,
    nug.min = 1e-8,
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
      #self$theta_update()
      self$all_update()
    },
    update_params = function () {#browser()
      self$K <- corr_gauss_mat(self$X, theta=self$theta) + diag(self$nug, self$N)
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
      kxx <- corr_gauss_mat(XX, theta=self$theta)
      kx.xx <- corr_gauss_mat(XX, self$X, theta=self$theta)

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
      K <- corr_gauss_mat(self$X, theta=theta) + diag(self$nug, self$N)
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
    },




    deviance2 = function (theta, nug) { # joint deviance
      K <- corr_gauss_mat(self$X, theta=theta) + diag(nug, self$N)
      Kchol <- try(chol(K))
      if (inherits(Kchol, "try-error")) {return(Inf)}
      Kinv <- chol2inv(Kchol)
      mu_hat <- sum(Kinv %*% self$Z) / sum(Kinv)
      s2_hat <- c(t(self$Z - self$mu_hat) %*% Kinv %*% (self$Z - mu_hat) / self$N)
      detK <- prod(diag(Kchol))^2
      log(detK) + self$N * log((self$Z - mu_hat) %*% solve(K, self$Z - mu_hat))
    },
    deviance_log2 = function (params) { # joint deviance
      beta <- params[-length(params)]
      nug <- params[length(params)]
      theta <- 10^beta
      self$deviance2(theta=theta, nug=nug)
    },
    deviance_search2 = function () {
      # Joint MLE deviance search with genetic alg, slow and bad
      #optimize(self$deviance_log, interval = c(-10,10), maximum = F, tol=1e-3)
      nvars <- self$D
      suppressWarnings(
        rgenoud::genoud(self$deviance_log2, nvars=nvars + 1,#self$D,#1,
                        Domains = rbind(matrix(c(-10,10), nvars, 2, byrow=T),c(0,1)), max = F,
                        solution.tolerance=1e-3,
                        pop.size=min(1000, 200*self$D),
                        max.generations = min(10, 3*self$D),
                        print.level = 0
        )
      )
    },
    deviance_search3 = function (restarts = 5) {#browser()
      # Joint MLE search with L-BFGS-B, with restarts
      #optimize(self$deviance_log, interval = c(-10,10), maximum = F, tol=1e-3)
      #nvars <- self$D
      lower <- c(rep(-10, self$D), self$nug.min)
      upper <- c(rep(10, self$D), Inf)
      # maybe start by calculating deviance for current params
      best <- try(
        optim(c(log(self$theta, 10),self$nug), self$deviance_log2, method="L-BFGS-B", lower=lower, upper=upper, hessian=F)
      )
      if (inherits(best, "try-error")) {
        best <- list(par=c(log(self$theta, 10), self$nug), value = self$deviance_log2(c(log(self$theta, 10))))
      }
      if (restarts >= 1) {
        for (i in 1:restarts) {
          current <- try(
            optim(c(log(self$theta, 10) + rnorm(self$D,0,2),self$nug + rexp(1, 1e4)), self$deviance_log2, method="L-BFGS-B", lower=lower, upper=upper, hessian=F)
          )
          if (!inherits(current, "try-error")) {
            if (current$value < best$value) {
              best <- current
            }
          }
        }
      }
      best
    },
    all_update = function () { # update theta and nugget
      pars <- self$deviance_search3()$par
      #self$theta <- 10 ^ self$deviance_search()$minimum
      self$nug <- pars[length(pars)]
      self$theta <- 10 ^ pars[-length(pars)]
      self$update_params()
    },

    deviance_searchnug = function() {
      optim(self$nug, function(nnug) {self$deviance2(theta=self$theta, nug=nnug)}, method="L-BFGS-B", lower=0, upper=Inf, hessian=F)$par
    },
    nugget_update = function () {
      nug <- self$deviance_searchnug()
      #self$theta <- 10 ^ self$deviance_search()$minimum
      self$nug <- nug
      self$update_params()
    },
    mean_grad = function (XX) {#browser()

      if (!is.matrix(XX)) {
        if (self$D == 1) XX <- matrix(XX, ncol=1)
        else if (length(XX) == self$D) XX <- matrix(XX, nrow=1)
        else stop('Predict input should be matrix')
      } else {
        if (ncol(XX) != self$D) {stop("Wrong dimension input")}
      }
      kx.xx <- t(corr_gauss_mat(XX, self$X, theta=self$theta))

      #grad <- numeric(self$D)
      #drdx <- matrix(NA, self$N, self$D)
      #for (i in 1:self$N) {
      #  for (j in 1:self$D) {
      #    drdx[i, j] <- -2 * self$theta[j] * (XX[j, 1] - self$X[i, j]) * kx.xx[i, 1]
      #  }
      #}
      #drdx <- -2 * outer(1:self$N, 1:self$D, Vectorize(function(i,j) {self$theta[j] * (XX[j, 1] - self$X[i, j]) * kx.xx[i, 1]}))
      #t(drdx) %*% self$Kinv %*% (self$Z - self$mu_hat)

      # mult points and dims
      #drdx <- vapply(1:nrow(XX), Vectorize(function(k) -2 * outer(1:self$N, 1:self$D, Vectorize(function(i,j) {self$theta[j] * (XX[k, j] - self$X[i, j]) * kx.xx[i, k]}))), numeric(self$N))
      #t(drdx) %*% self$Kinv %*% (self$Z - self$mu_hat)

      grad <-   vapply(1:nrow(XX),

                      Vectorize(
                        function(k) {
                          t(-2 * outer(1:self$N, 1:self$D, Vectorize(function(i,j) {self$theta[j] * (XX[k, j] - self$X[i, j]) * kx.xx[i, k]}))
                          )  %*%self$Kinv %*% (self$Z - self$mu_hat)
                        }
                      )

               #, matrix(0, self$N, nrow(XX)) #rep(0, self$N)
               , numeric(self$D) #rep(0, self$N)

        )
      if (self$D == 1) return(grad)
      t(grad)

    },
    mean_grad_norm = function (XX) {#browser()
      grad <- self$mean_grad(XX)
      if (!is.matrix(grad)) return(abs(grad))
      apply(grad,1, function(xx) {sqrt(sum(xx^2))})
    }
  ),
  private = list(

  )
)
