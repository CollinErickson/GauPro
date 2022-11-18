library(testthat)

context("Test kernel model LOO")

test_that("kernel model LOO works", {
  # Check if LOO predictions match actual on banana function,  i.e. check shortcut
  set.seed(0)
  n <- 80
  d <- 2
  f1 <- function(x) {abs(sin(2*pi*x[1])) + x[2]^2}
  X1 <- matrix(runif(n*d),n,d)
  Z1 <- apply(X1,1,f1) * 9.3  + rnorm(n, 0, 1e0)
  expect_error(gp <- GauPro_kernel_model_LOO$new(X=X1, Z=Z1, kernel="matern52",
                                                 nug.min=1e-6),
               NA)
  # ContourFunctions::cf(gp$predict, pts=X1)
  nn <- 1e3
  XX <- matrix(runif(nn*d),nn,d)
  ZZ <- apply(XX, 1, f1) * 9.3
  gp$use_LOO <- T
  # Predict
  expect_error(ZZhat <- gp$predict(XX, se=F), NA)
  expect_error(ZZhat <- gp$predict(XX, se=T), NA)
  # Predict mean
  expect_error(ZZhat <- gp$predict(XX, se=T, mean_dist=T), NA)
})

test_that("pred_LOO", {
  # Test pred_LOO
  n <- 80
  d <- 2
  f1 <- function(x) {abs(sin(2*pi*x[1])) + x[2]^2}
  X1 <- matrix(runif(n*d),n,d)
  Z1 <- apply(X1,1,f1) * 9.3 + rnorm(n, 0, 1e0)
  expect_error(gp <- GauPro_kernel_model$new(X=X1, Z=Z1, kernel="matern52",
                                             nug.min=1e-6),
               NA)

  # Make sure that the LOO predictions match LOO
  expect_error(ZLOO <- gp$pred_LOO(se=T), NA)
  gp2 <- gp$clone(deep=T)
  loo_means <- numeric(n)
  loo_ses <- numeric(n)
  for (i in 1:n) {
    gpi <- gp$clone(deep=T);
    expect_error(gpi$update(Xall=X1[-i,],Zall=Z1[-i], no_update = TRUE), NA)
    # if (T) { #set mu and s2 back to original values
    #   # This makes differences ~ 1e-15 instead of 1e-4, not sure if it is recommended though
    #   gpi$s2_hat <- gp$s2_hat
    #   gpi$mu_hat <- gp$mu_hat
    # }
    expect_error(gpp <- gpi$predict(X1[i,],se=T), NA)
    loo_means[i] <- gpp$me
    loo_ses[i] <- gpp$se
  }
  # cbind(ZLOO$fit, loo_means)
  # summary(ZLOO$fit - loo_means)
  # plot(ZLOO$fit, loo_means)

  expect_true(max(abs(ZLOO$fit - loo_means)) < 1e-4)
})
