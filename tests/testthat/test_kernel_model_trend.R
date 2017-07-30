test_that("trend_0 works", {
  set.seed(0)
  n <- 20
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
  y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Matern32$new(1), trend=trend_0$new(), parallel=FALSE, verbose=10, nug.est=T, nug.min=0)

  expect_equal(gp$kernel$beta, 0.4609827, tolerance=.01)
  expect_equal(gp$kernel$s2, 1.083443, tolerance=.01)
  expect_equal(log(gp$nug,10), -9.710726, tolerance=.1)
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:2], nuglog=x[3])}, x=c(-.7,.227, -3.66)),
    c(-877.52290, -643.49553, -30.76981),
    gp$deviance_grad(params = c(-.7,.227), nug.update=T, nuglog=-3.66, trend_update=TRUE),
    tol=.001
  )
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:2], nuglog=x[3])}, x=c(2.5,-.2, -5)),
    c(2.393571e+01, 3.066583e+01, 7.917006e-04),
    gp$deviance_grad(params = c(2.5,-.2), nug.update=T, nuglog = -5, trend_update=TRUE),
    tol=.01
  )
})

test_that("trend_c works", {
  set.seed(0)
  n <- 20
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
  y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Matern32$new(1), trend=trend_c$new(), parallel=FALSE, verbose=10, nug.est=T, nug.min=0)

  expect_equal(gp$trend$m, -0.002844272, tolerance=.01)
  expect_equal(gp$kernel$beta, 0.4609827, tolerance=.01)
  expect_equal(gp$kernel$s2, 1.083443, tolerance=.01)
  expect_equal(log(gp$nug,10), -11.37841, tolerance=.1)
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[2:3], nuglog=x[4],trend_params=x[1])}, x=c(10.2,-.7,.227, -3.66)),
    c(14.26316, -864.14491, -811.75501,  -32.01259),
    gp$deviance_grad(params = c(-.7,.227), nug.update=T, nuglog=-3.66, trend_params=10.2, trend_update=TRUE),
    tol=.001
  )
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[2:3], nuglog=x[4],trend_params=x[1])}, x=c(-4.5, 2.5,-.2, -5)),
    c(-1.221971e+02,  2.985039e+02, -6.024080e+02, -2.102204e-03),
    gp$deviance_grad(params = c(2.5,-.2), nug.update=T, nuglog = -5, trend_params=-4.5, trend_update=TRUE),
    tol=.01
  )
})

test_that("trend_LM works", {
  set.seed(0)
  n <- 20
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
  y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Matern32$new(1), trend=trend_LM$new(D=1), parallel=FALSE, verbose=10, nug.est=T)

  expect_equal(gp$trend$b, -0.5470029, tolerance=.01)
  expect_equal(gp$trend$m, 1.087638, tolerance=.01)
  expect_equal(gp$kernel$beta, 0.4007064, tolerance=.01)
  expect_equal(gp$kernel$s2, 1.247229, tolerance=.01)
  expect_equal(log(gp$nug,10), -13.01126, tolerance=1)
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[3:4], nuglog=x[5],trend_params=x[1:2])}, x=c(10.2, 2,-.7,.227, -3.66)),
    c(15.655129,   2.753122, -822.992291, -806.199697,  -29.333112),
    gp$deviance_grad(params = c(-.7,.227), nug.update=T, nuglog=-3.66, trend_params=c(10.2,2), trend_update=TRUE),
    tol=.001
  )
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[3:4], nuglog=x[5],trend_params=x[1:2])}, x=c(-4.5, -3.2, 2.5,-.2, -5)),
    c(-1.656454e+02, -8.804362e+01,  5.253556e+02, -1.137605e+03, -4.612252e-03),
    gp$deviance_grad(params = c(2.5,-.2), nug.update=T, nuglog = -5, trend_params=c(-4.5,-3.2), trend_update=TRUE),
    tol=.01
  )
})
