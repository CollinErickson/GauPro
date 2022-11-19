library(testthat)

printkern <- interactive()

test_that("All kernels work in 1-D", {
  kern_chars <- c("gauss", "exp", "m32", "m52",
                  "ratquad", "periodic", "cubic",
                  "triangle", "white")
  kernels <- list(Gaussian$new(0),
                  Exponential$new(0),
                  Matern32$new(0),
                  Matern52$new(0),
                  RatQuad$new(0,0),
                  Periodic$new(p=1,alpha=1),
                  Cubic$new(0),
                  Triangle$new(0),
                  White$new(s2=.65))
  expect_equal(length(kern_chars), length(kernels))
  set.seed(0)
  n <- 20
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
  y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  x1 <- .34343
  x2 <- matrix(c(.2352,.4745,.34625,.97654,.16435,.457, .354, .976,.234, .623), ncol=1)
  for (ikern in seq_along(kernels)) {
    kern_char <- kern_chars[ikern]
    kernel <- kernels[[ikern]]

    if (exists('printkern') && printkern) cat("1D", j, kern_char, "\n")
    set.seed(0)

    # Fit GP using kernel
    #' kernel=RatQuad$new(0.1,0.1)
    expect_error({
      gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=kernel,
                                    parallel=FALSE, verbose=0, nug.est=T,
                                    restarts=0)
    }, NA)

    x1_C <- gp$kernel$k(x=x1)
    expect_is(object = x1_C, class = 'numeric', info = class(kernel)[1])
    expect_equal(object = length(x1_C), expected = 1, info = class(kernel)[1])

    x2_C <- gp$kernel$k(x=x2)
    expect_is(object = x2_C, class = 'matrix', info = class(kernel)[1])
    # expect_length(object = x2_C, n = 100, info = class(kernel)[1])
    expect_equal(object = length(x2_C), expected = 100, info = class(kernel)[1])

    # Check C_dC
    C <- gp$kernel$k(x=x)
    dC <- gp$kernel$dC_dparams(X=x, nug=gp$nug)
    C_dC <- gp$kernel$C_dC_dparams(X=x, nug=gp$nug)
    expect_equal(object = C+diag(gp$s2_hat * gp$nug, n), expected = C_dC[[1]], info = class(kernel)[1])
    expect_equal(object = dC, expected = C_dC[[2]], info = class(kernel)[1])
  }
})

test_that("All kernels work in 2-D", {
  kern_chars <- c("gauss", "exp", "m32", "m52",
                  "ratquad", "cubic", "triangle", "periodic")
  kernels <- list(Gaussian$new(c(0,.3)),
                  Exponential$new(c(0,.3)),
                  Matern32$new(c(0,.3)),
                  Matern52$new(c(0,.3)),
                  RatQuad$new(c(0,.3),0),
                  Cubic$new(D=2),
                  Triangle$new(D=2),
                  Periodic$new(p=c(.4,.3),alpha=1))
  expect_equal(length(kern_chars), length(kernels))
  set.seed(0)
  n <- 30
  x <- matrix(runif(2*n), ncol=2)
  f <- function(x) {sin(2*pi*x[1]) + .001*sin(8*pi*x[1]) + exp(x[2]) + .03*x[1]*x[2] +rnorm(1,0,.03)}
  y <- apply(x, 1, f) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  x1 <- c(.34343, .65)
  x2 <- matrix(c(.2352,.4745,.34625,.97654,.16435,.457, .354, .976,.234, .623), ncol=2)
  for (ikern in seq_along(kernels)) {
    kern_char <- kern_chars[ikern]
    kernel <- kernels[[ikern]]

    if (exists('printkern') && printkern) cat("2D", j, kern_char, "\n")
    set.seed(0)

    # Fit GP using kernel
    #' kernel=RatQuad$new(0.1,0.1)
    expect_error({
      gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=kernel,
                                    parallel=FALSE, verbose=0, nug.est=T,
                                    restarts=0)
    }, NA)

    x1_C <- gp$kernel$k(x=x1)
    expect_is(object = x1_C, class = 'numeric', info = class(kernel)[1])
    expect_equal(object = length(x1_C), expected = 1, info = class(kernel)[1])

    x2_C <- gp$kernel$k(x=x2)
    expect_is(object = x2_C, class = 'matrix', info = class(kernel)[1])
    expect_equal(object = length(x2_C), expected = 25, info = class(kernel)[1])

    # Check C_dC
    C <- gp$kernel$k(x=x)
    dC <- gp$kernel$dC_dparams(X=x, nug=gp$nug)
    C_dC <- gp$kernel$C_dC_dparams(X=x, nug=gp$nug)
    expect_equal(object = C+diag(gp$s2_hat * gp$nug, n), expected = C_dC[[1]], info = class(kernel)[1])
    expect_equal(object = dC, expected = C_dC[[2]], info = class(kernel)[1])
  }
})

test_that("Gaussian kernel has correct grad", {
  n <- 30
  d <- 3
  X <- matrix(runif(d*n), ncol=d)
  # f <- function(x) {sin(2*pi*x[1]) + .001*sin(8*pi*x[1]) + exp(x[2]) + .03*x[1]*x[2] +rnorm(1,0,.03)}
  # y <- apply(x, 1, f) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  # x1 <- c(.34343, .65)
  # x2 <- matrix(c(.2352,.4745,.34625,.97654,.16435,.457, .354, .976,.234, .623), ncol=2)
  k1 <- Gaussian$new(D=d, beta=runif(d,-1,1), s2=rexp(1))
  # gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=kernel, parallel=FALSE, verbose=0, nug.est=T, restarts=1)
  eps <- 1e-8
  # Check correlation parameters
  for (idim in 1:d) {
    epsvec <- rep(0, d)
    epsvec[idim] <- eps
    # k1 <- kernel#$k(x)
    CdC0 <- k1$C_dC_dparams(params = c(k1$beta-epsvec/2, k1$logs2), X=X, nug=1e-4)
    CdC1 <- k1$C_dC_dparams(params = c(k1$beta+epsvec/2, k1$logs2), X=X, nug=1e-4)
    d1 <- (CdC1$C - CdC0$C) / eps
    d2 <- CdC0$dC_dparams[idim,,]
    d1
    d2
    expect_equal(d1, d2, tolerance = 1e-6)
  }
  # Check logs2
  CdC0 <- k1$C_dC_dparams(params = c(k1$beta, k1$logs2-eps/2), X=X, nug=1e-4)
  CdC1 <- k1$C_dC_dparams(params = c(k1$beta, k1$logs2+eps/2), X=X, nug=1e-4)
  d1 <- (CdC1$C - CdC0$C) / eps
  d2 <- CdC0$dC_dparams[d+1,,]
  d1
  d2
  expect_equal(d1, d2, tolerance = 1e-6)
})

test_that("Cubic kernel has correct grad", {
  n <- 30
  d <- 3
  X <- matrix(runif(d*n), ncol=d)
  k1 <- Cubic$new(D=d, beta=runif(d,-1,1), s2=rexp(1))
  eps <- 1e-8
  # Check correlation parameters
  for (idim in 1:d) {
    epsvec <- rep(0, d)
    epsvec[idim] <- eps
    # k1 <- kernel#$k(x)
    CdC0 <- k1$C_dC_dparams(params = c(k1$beta-epsvec/2, k1$logs2), X=X, nug=1e-4)
    CdC1 <- k1$C_dC_dparams(params = c(k1$beta+epsvec/2, k1$logs2), X=X, nug=1e-4)
    d1 <- (CdC1$C - CdC0$C) / eps
    d2 <- CdC0$dC_dparams[idim,,]
    d1
    d2
    expect_equal(d1, d2, tolerance = 1e-6)
  }
  # Check logs2
  CdC0 <- k1$C_dC_dparams(params = c(k1$beta, k1$logs2-eps/2), X=X, nug=1e-4)
  CdC1 <- k1$C_dC_dparams(params = c(k1$beta, k1$logs2+eps/2), X=X, nug=1e-4)
  d1 <- (CdC1$C - CdC0$C) / eps
  d2 <- CdC0$dC_dparams[d+1,,]
  d1
  d2
  expect_equal(d1, d2, tolerance = 1e-6)
})



test_that("FactorKernel has correct grad", {
  kfac <- FactorKernel$new(D=1, nlevels = 3, xindex = 1, s2 = .333)
  kfac$p <- c(.1,.2,.3)
  X <- matrix(c(1,2,3,1,2,3), ncol=1)
  eps <- 1e-8
  for (i in 1:length(kfac$p)) {
    # print(i)
    epsvec <- rep(0, length(kfac$p))
    epsvec[i] <- eps
    CdC0 <- kfac$C_dC_dparams(params = c(kfac$p, kfac$logs2), X=X, nug=1e-4)
    CdC1 <- kfac$C_dC_dparams(params = c(kfac$p+epsvec, kfac$logs2), X=X, nug=1e-4)
    d1 <- (CdC1$C - CdC0$C) / (1e-8)
    d2 <- CdC0$dC_dparams[i,,]
    d1
    d2
    expect_equal(d1, d2)
  }
})



test_that("OrderedFactorKernel has correct grad", {
  kfac <- OrderedFactorKernel$new(D=1, nlevels = 3, xindex = 1, s2 = 1.333)
  kfac$p <- c(1.1,1.2)
  X <- matrix(c(1,2,3,1,2,3), ncol=1)
  eps <- 1e-6
  for (i in 1:length(kfac$p)) {
    # print(i)
    epsvec <- rep(0, length(kfac$p))
    epsvec[i] <- eps
    CdC0 <- kfac$C_dC_dparams(params = c(kfac$p-epsvec, kfac$logs2), X=X, nug=1e-4)
    CdC1 <- kfac$C_dC_dparams(params = c(kfac$p+epsvec, kfac$logs2), X=X, nug=1e-4)
    d1 <- (CdC1$C - CdC0$C) / (2*eps)
    d2 <- CdC0$dC_dparams[i,,]
    d1
    d2
    d1/d2
    expect_equal(d1, d2, tolerance = 1e-4)
  }
})

# Doesn't work since find_kfd isn't exported
if (F) {
  test_that('ignore inds', {
    k <- IgnoreIndsKernel$new(ignoreinds = c(1,3,5),
                              k = Gaussian$new(D=2) *
                                OrderedFactorKernel$new(D=2, xindex = 2, nlevels = 4)) *
      FactorKernel$new(nlevels=11, xindex=1, D=5) *
      OrderedFactorKernel$new(nlevels=3, xindex=3, D=5) *
      LatentFactorKernel$new(nlevels=5, xindex=5, D=5)
    k
    fk <- find_kernel_factor_dims(k)
    expect_equal(fk, c(4,4,1,11,3,3,5,5))
  })
}
