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

  expect_equal(gp$kernel$beta, 1.804936, tolerance=.1)
  expect_equal(gp$kernel$s2, .895724, tolerance=.1)
  expect_equal(gp$nug, 0.0000020596, tolerance=.0001)
  expect_equal(
    # numDeriv::grad(func = gp$deviance, x=c(5,1)),
    c(6.799652, 11.370354),
    gp$deviance_grad(params = c(3,1), nug.update=F),
    tol=.01
  )
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:2], nuglog=x[3])}, x=c(2.5,.7, -5)),
    c(30.00001, 17.82337, 6.405043e-4),
    gp$deviance_grad(params = c(2.5,.7), nug.update=T, nuglog = -5),
    tol=.01
  )
})
