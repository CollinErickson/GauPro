test_that("All kernels work in 1-D", {
  kernels <- list(Gaussian$new(0),
                  Exponential$new(0),
                  Matern32$new(0),
                  Matern52$new(0),
                  RatQuad$new(0,0),
                  Periodic$new(p=1,alpha=1),
                  White$new(s2=.65))
  set.seed(0)
  n <- 20
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
  y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  x1 <- .34343
  x2 <- matrix(c(.2352,.4745,.34625,.97654,.16435,.457, .354, .976,.234, .623), ncol=1)
  for (kernel in kernels) {
    set.seed(0)

    # Fit GP using kernel
    #' kernel=RatQuad$new(0.1,0.1)
    gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=kernel, parallel=FALSE, verbose=0, nug.est=T, restarts=1)

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
    #' kernel=RatQuad$new(0.1,0.1)
    gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=kernel, parallel=FALSE, verbose=0, nug.est=T, restarts=1)

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
