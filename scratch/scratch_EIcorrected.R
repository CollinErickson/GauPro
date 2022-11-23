# EI with high noise. Gives too much improvement.
d <- 1
n <- 20
# x <- runif(n)
x <- seq(0,1,l=n)
x <- c(runif(n/2,0,.3), runif(5,.3,.7), runif(n/2-5,.7,1))
y <- sin(2*pi*x) + rnorm(length(x), 0, .2)
plot(x, y)
# gp1 <- GauPro_kernel_model$new(y ~ x)
gp1 <- GauPro_kernel_model$new(x, y, kernel='m52')
gp1$cool1Dplot()
gp1$plot1D()
curve(gp1$pred(x, se.fit = T)$se, ylim=c(0,.3))
curve(gp1$pred(x, se.fit = T, mean_dist = T)$se, add=T, col=2)
# EI is way too big
gp1$maxEI(minimize = T)
curve(sapply(x, function(x)gp1$EI(x, minimize = T)))


# What should it be?
# If you add the point.
# Sample over dist.
# Add that sample to the model.
# Find new mean prediction at that point.
# KG runs optimization after that.


EI2 <- function(self, x, minimize=FALSE) {
  self <- gp1
  predX <- self$pred(self$X)

  x <- .35
  X <- self$X
  Z <- self$Z
  u <- x
  a <- x

  mu_a <- self$trend$Z(a)
  mu_X <- c(self$trend$Z(X))
  KX <- self$K
  KXinv <- self$Kinv
  KuX <- matrix(self$kernel$k(u, X), nrow=1)
  KXu <- t(KuX)
  KaX <- self$kernel$k(a, X)
  KXa <- t(KaX)
  Kau <- self$kernel$k(a, u)
  Ka <- self$kernel$k(a)
  Ku <- self$kernel$k(u) + self$s2_hat*self$nug
  Ma <- c(KaX, Kau)
  Mcinv <- rbind(cbind(KX, KXu),
                 cbind(KuX, Ku))
  Mc <- solve(Mcinv)
  KugivenX <- Ku - KuX %*% KXinv %*% KXu

  self$pred(a, se=T, mean_dist = T)

  mu_m <- c(mu_a + Ma %*% Mc %*% rbind(Z-mu_X, c(KuX%*%KXinv%*%(Z-mu_X))))
  mu_m

  cov_m <- matrix(Ma, nrow=1) %*% Mc %*%
    rbind(cbind(matrix(0, ncol=ncol(KX)+1, nrow=nrow(KX))),
                               cbind(matrix(0, ncol=ncol(KX),   nrow=1),
                                     KugivenX)) %*%
    Mc %*% matrix(Ma, ncol=1)
  cov_m
  # This is right! Matches simulation below.


  # Sample to estimate dist
  ns <- 1e3
  ms <- rep(0, ns)
  for (i in 1:ns) {
    gpi <- self$clone(deep=T)
    zi <- gpi$sample(u)
    gpi$update(Xnew=u, Znew=zi, no_update = T)
    ms[i] <- gpi$pred(u)
  }
  summary(ms)
  var(ms)
  hist(ms, freq=F)
  curve(dnorm(x, mu_m, sqrt(c(cov_m))), add=T, col=2)
}
