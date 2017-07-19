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

  expect_equal(gp$kernel$theta, 4.826752, tolerance=.1)
  expect_equal(gp$kernel$s2, 1.963462, tolerance=.1)
  expect_equal(gp$nug, 0.000435802, tolerance=.01)
  expect_equal(
    # numDeriv::grad(func = function(x) gp$deviance(params=x[1:2], nuglog=x[3]), x=c(5,1, -4)),
    c(-6.641381,  -98.038047, -213.498400),
    gp$deviance_grad(params = c(5,1), nug.update=T, nuglog=-4),
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

  expect_equal(gp$kernel$beta, 1.185298, tolerance=.1)
  expect_equal(gp$kernel$s2, 0.3809001, tolerance=.1)
  expect_equal(gp$nug, 0.001037888, tolerance=.0001)
  expect_equal(
    # numDeriv::grad(func = function(x) gp$deviance(params=x[1:2], nuglog=x[3]), x=c(.2,1.1, -4.5)),
    c(-272.12230, -102.95356,  -63.01746),
    gp$deviance_grad(params = c(.2,1.1), nug.update=T, nuglog=-4.5),
    tol=.001
  )
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:2], nuglog=x[3])}, x=c(2.5,-.2, -5)),
    c(30.89216, 26.80619,  0.0006299453),
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

  expect_equal(gp$kernel$l, 0.1806493, tolerance=.01)
  expect_equal(gp$kernel$s2, 0.3809001, tolerance=.01)
  expect_equal(gp$nug, 0.001037888, tolerance=.0001)
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:2], nuglog=x[3])}, x=c(-.7,.227, -3.66)),
    c(-20.509641,  11.177885,  -1.592343),
    gp$deviance_grad(params = c(-.7,.227), nug.update=T, nuglog=-3.66),
    tol=.001
  )
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:2], nuglog=x[3])}, x=c(2.5,-.2, -5)),
    c(3847404.2, -2317793.3, -1894093.7),
    gp$deviance_grad(params = c(2.5,-.2), nug.update=T, nuglog = -5),
    tol=100
  )
})

