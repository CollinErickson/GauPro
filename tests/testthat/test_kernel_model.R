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
