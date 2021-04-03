test_that("1D data works", {

  n <- 12
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  y <- sin(2*pi*x) + rnorm(n,0,1e-1)
  gp <- GauPro(X=x, Z=y, parallel=FALSE)
  expect_that(gp, is_a("GauPro"))
  expect_that(gp, is_a("R6"))
})


test_that("corr works", {
  m1 <- outer(1:10, 1:10, Vectorize(function(i,j) {exp(-sum((1e-2) * (i-j-5)^2))}))
  m2 <- corr_gauss_matrixC(matrix(1:10,ncol=1), matrix(6:15,ncol=1), 1e-2)
  m3 <- corr_gauss_matrix(matrix(1:10,ncol=1), matrix(6:15,ncol=1), 1e-2)
  expect_equal(m1, m2)
  expect_equal(m1, m3)
  expect_equal(m2, m3)
})
