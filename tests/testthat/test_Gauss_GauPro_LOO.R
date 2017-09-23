context("Test GauPro_Gauss_LOO")

test_that("GauPro_Gauss_LOO works", {
  # Check if LOO predictions match actual on banana function,  i.e. check shortcut
  set.seed(0)
  n <- 80
  d <- 2
  f1 <- function(x) {abs(sin(2*pi*x[1])) + x[2]^2}
  X1 <- matrix(runif(n*d),n,d)
  Z1 <- apply(X1,1,f1) * 9.3 # + rnorm(n, 0, 1e-3)
  gp <- GauPro_Gauss_LOO$new(X=X1, Z=Z1)
  # ContourFunctions::cf(gp$predict, pts=X1)
  nn <- 1e3
  XX <- matrix(runif(nn*d),nn,d)
  ZZ <- apply(XX, 1, f1) * 9.3
  gp$use_LOO <- T
  ZZhat <- gp$predict(XX, se=T)

  ZLOO <- gp$pred_LOO(se=T)
  gp2 <- gp$clone(deep=T)
  loo_means <- numeric(n)
  loo_ses <- numeric(n)
  for (i in 1:n) {
    gpi <- gp$clone(deep=T);
    gpi$update(Xall=X1[-i,],Zall=Z1[-i], no_update = TRUE);
    if (T) { #set mu and s2 back to original values
      # This makes differences ~ 1e-15 instead of 1e-4, not sure if it is recommended though
      gpi$s2_hat <- gp$s2_hat
      gpi$mu_hat <- gp$mu_hat
    }
    gpp <- gpi$predict(X1[i,],se=T)
    loo_means[i] <- gpp$me
    loo_ses[i] <- gpp$se
  }
  # cbind(ZLOO$fit, loo_means)
  # summary(ZLOO$fit - loo_means)

  expect_true(max(abs(ZLOO$fit - loo_means)) < 1e-8)
})
