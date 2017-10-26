# Test to check that grad returns right length and
# agrees with numerical grad.
test_that("kernel grad works", {




  kernels <- list(Gaussian$new(c(0,.3)),
                  Exponential$new(c(0,.3)),
                  Matern32$new(c(0,.3)),
                  Matern52$new(c(0,.3)),
                  RatQuad$new(c(0,.3),0),
                  Periodic$new(p=c(.4,.3),alpha=1),
                  White$new(s2=.96))
  set.seed(0)
  n <- 30
  x <- matrix(runif(2*n), ncol=2)
  f <- function(x) {sin(2*pi*x[1]) + .001*sin(8*pi*x[1]) + exp(x[2]) + .03*x[1]*x[2] +rnorm(1,0,.03)}
  y <- apply(x, 1, f) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  x1 <- c(.34343, .65)
  x2 <- matrix(c(.2352,.4745,.34625,.97654,.16435,.457, .354, .976,.234, .623), ncol=2)
  for (kernel in kernels) {
    set.seed(0)

    # Fit GP using kernel
    #' kernel=Gaussian$new(c(0.1,.2))
    gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=kernel, parallel=FALSE, verbose=10, nug.est=T, restarts=0)

    grad1 <- gp$grad(x1)
    expect_is(object = grad1, class = 'matrix')
    expect_length(object = grad1, n = 2)

    grad2 <- gp$grad(x2)
    expect_is(object = grad2, class = 'matrix')
    expect_length(object = grad2, n = 10)

    # Check grad with numerical grad
    if (requireNamespace("numDeriv", quietly = TRUE)) {
      expect_equal(c(gp$grad(x1)), c(numDeriv::grad(gp$predict, x1)), tol=.1)
    }
  }
})

test_that("grad_norm2 for Gaussian", { # only implemented for Gaussian now
  set.seed(0)
  n <- 20
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  f <- Vectorize(function(x) {sin(2*pi*x) + .5*sin(4*pi*x) +rnorm(1,0,.3)}+10*x)
  y <- 123 + f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian$new(beta=0.152, s2=10^1.7194, beta_est=F, s2_est=F), restarts=0, no_update=T, parallel=FALSE, verbose=10, nug.est=F, param.est=F, nug=0.00357)
  expect_equal(unlist(gp$grad_norm2_dist(matrix(.1,ncol=1))), c(mean=324.1904, var=3781.708), tol=1)
  set.seed(1)
  ts <- gp$grad_norm2_sample(matrix(.1,ncol=1), n=1e4)
  expect_equal(dim(ts), c(1, 1e4))
  expect_equal(c(mean(ts), sd(ts)^2), c(325.7848, 3877.7587), tol=1)
  gd <- gp$grad_dist(XX=matrix(.1))
  expect_is(gd, "list")
  expect_length(gd, 2)
  expect_equal(gd[[1]], array(data = 17.9775, dim = c(1,1)), tol=.1)
  expect_equal(gd[[2]], array(data = 2.923747, dim = c(1,1,1)), tol=.1)

  # Check 2D
  set.seed(0)
  n <- 30
  x <- matrix(runif(2*n), ncol=2)
  f <- function(x) {sin(2*pi*x[1]) + .5*sin(4*pi*x[1]) +rnorm(1,0,.03) + x[2]^2}
  y <- apply(x, 1, f) #f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian$new(beta=c(.966, -.64),s2=10^.4848), parallel=FALSE, verbose=10, nug.est=F, nug=0.0001721, param.est=F)
  # Check grad_norm2_dist is correct
  tx <- matrix(c(.1,.2),ncol=2)
  expect_equal(unlist(gp$grad_norm2_dist(tx)), c(mean=45.46676, var=31.37130), tol=.0001)

  # Check that it matches samples
  set.seed(1)
  ts <- gp$grad_norm2_sample(tx, n=1e4)
  expect_equal(dim(ts), c(1, 1e4))
  # expect_equal(c(mean(ts), sd(ts)^2), c(47.46267, 19.65818), tol=.001)
  expect_equal(unlist(gp$grad_norm2_dist(tx)), c(mean=mean(ts), var=var(c(ts))), tol=.01)

  # Alternate way to calculate just the mean
  tmp2 <- sum(gp$grad(tx)^2) +sum(eigen(gp$grad_dist(tx)$cov[1,,])$val)
  expect_equal(gp$grad_norm2_dist(tx)$mean, tmp2, tol=.00001)

  gd <- gp$grad_dist(XX=matrix(c(.1,.2), ncol=2))
  expect_is(gd, "list")
  expect_length(gd, 2)
  expect_equal(dim(gd[[1]]), c(1,2))
  expect_equal(dim(gd[[2]]), c(1,2,2))

})
