test_that("1D data works", {

  n <- 12
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  y <- sin(2*pi*x) + rnorm(n,0,1e-1)
  gp <- GauPro(X=x, Z=y, parallel=FALSE)
  expect_that(gp, is_a("GauPro"))
  expect_that(gp, is_a("R6"))
})
