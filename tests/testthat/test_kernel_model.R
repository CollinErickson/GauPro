test_that("kernel_model works", {

  set.seed(0)
  n <- 20
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
  y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian$new(1), parallel=FALSE, verbose=10, nug.est=T)
  # gp$cool1Dplot()
  # numDeriv::grad(func = gp$deviance, x=c(5,1))
  # gp$deviance_grad(params = c(5,1), nug.update=F)

  expect_equal(gp$kernel$theta, 4.317251, tolerance=.1)
  expect_equal(gp$kernel$s2, 2.223157, tolerance=.1)
  expect_equal(gp$nug, 0.00085965, tolerance=.01)
  expect_equal(
    # numDeriv::grad(func = gp$deviance, x=c(5,1)),
    c(-1.079865, -2.202628),
    gp$deviance_grad(params = c(5,1), nug.update=F),
    tol=.01
    )
})
test_that("kernel_Gaussian_beta works", {

  set.seed(0)
  n <- 20
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
  y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian_beta$new(1), parallel=FALSE, verbose=10, nug.est=T)
  # gp$cool1Dplot()
  # numDeriv::grad(func = gp$deviance, x=c(5,1))
  # gp$deviance_grad(params = c(5,1), nug.update=F)

  expect_equal(gp$kernel$beta, 0.8360693, tolerance=.1)
  expect_equal(gp$kernel$s2, 0.6292807, tolerance=.1)
  expect_equal(gp$nug, 0.000937203, tolerance=.0001)
  expect_equal(
    # numDeriv::grad(func = gp$deviance, x=c(5,1)),
    c(8.344209, 16.450576),
    gp$deviance_grad(params = c(1,0), nug.update=F),
    tol=.001
  )
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:2], nuglog=x[3])}, x=c(2.5,-.2, -5)),
    c(30.948458901, 26.684818099,  0.000629267),
    gp$deviance_grad(params = c(2.5,-.2), nug.update=T, nuglog = -5),
    tol=.001
  )
})
test_that("kernel_Gaussian_l works", {
  set.seed(0)
  n <- 20
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
  y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian_l$new(1), parallel=FALSE, verbose=10, nug.est=T)
  # gp$cool1Dplot()
  # numDeriv::grad(func = gp$deviance, x=c(5,1))
  # gp$deviance_grad(params = c(5,1), nug.update=F)

  expect_equal(gp$kernel$l, 0.2001134, tolerance=.01)
  expect_equal(gp$kernel$s2, 1.692982, tolerance=.01)
  expect_equal(gp$nug, 0.0002266279, tolerance=.0001)
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:2], nuglog=x[3])}, x=c(-.7,.227, -3.66)),
    c(-20.650727,  11.182670,  -1.593023),
    gp$deviance_grad(params = c(-.7,.227), nug.update=T, nuglog=-3.66),
    tol=.001
  )
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:2], nuglog=x[3])}, x=c(2.5,-.2, -5)),
    c(30.948458901, 26.684818099,  0.000629267),
    gp$deviance_grad(params = c(2.5,-.2), nug.update=T, nuglog = -5),
    tol=.001
  )
})

