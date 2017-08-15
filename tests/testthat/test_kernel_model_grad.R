# Test to check that grad returns right length and
# agrees with numerical grad.
test_that("kernel grad works", {




  kernels <- list(Gaussian$new(c(0,.3)),
                  Exponential$new(c(0,.3)),
                  Matern32$new(c(0,.3)),
                  Matern52$new(c(0,.3)),
                  RatQuad$new(c(0,.3),0),
                  Periodic$new(p=c(.4,.3),alpha=1))
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

    # Check C_dC

    expect_equal(c(gp$grad(x1)), c(numDeriv::grad(gp$predict, x1)), tol=.01)
  }
})
