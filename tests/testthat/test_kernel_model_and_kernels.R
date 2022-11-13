library(testthat)

set.seed(Sys.time())

printkern <- interactive()
# cat('printkern is', printkern, '\n')

# kernels work and have correct grads ----
test_that("kernels work and have correct grads", {
  n <- 20
  d <- 2
  x <- matrix(runif(n*d), ncol=d)
  f <- function(x) {abs(sin(x[1]^.8*6))^1.2 + log(1+x[2]) + x[1]*x[2]}
  y <- apply(x, 1, f) + rnorm(n,0,1e-4) #f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  kern_chars <- c('Gaussian', 'Matern32', 'Matern52',
                  'Triangle', 'Cubic', 'White',
                  'PowerExp', 'Periodic', "Exponential", "RatQuad",
                  "Ignore", "Product", "Sum")
  kern_list <- list(0,0,0,0,0,0,
                    0,0,0,0,
                    IgnoreIndsKernel$new(Gaussian$new(D=1), 2),
                    Gaussian$new(D=2)*PowerExp$new(D=2),
                    Cubic$new(D=2) * Triangle$new(D=2))
  stopifnot(length(kern_chars) == length(kern_list))
  for (j in 1:length(kern_chars)) {
    kern_char <- kern_chars[j]
    if (exists('printkern') && printkern) cat(j, kern_char, "\n")
    if (is.numeric(kern_list[[j]])) {
      kern <- eval(parse(text=kern_char))
      expect_is(kern, "R6ClassGenerator")
    } else {
      kern <- kern_list[[j]]
    }

    gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=kern, parallel=FALSE,
                                  verbose=0, nug.est=T, restarts=0)
    expect_is(gp, "GauPro")
    expect_is(gp, "R6")

    # Check predict
    pred1 <- predict(gp, runif(2))
    expect_true(is.numeric(pred1))
    expect_true(!is.matrix(pred1))
    expect_equal(length(pred1), 1)
    pred2 <- predict(gp, runif(2), se.fit=T)
    expect_true(is.data.frame(pred2))
    expect_equal(dim(pred2), c(1,3))
    expect_equal(colnames(pred2), c('mean', 's2', 'se'))
    pred3 <- predict(gp, matrix(runif(12), ncol=2))
    expect_true(is.numeric(pred3))
    expect_true(!is.matrix(pred3))
    expect_equal(length(pred3), 6)
    pred4 <- predict(gp, matrix(runif(12), ncol=2), se.fit=T)
    expect_true(is.data.frame(pred4))
    expect_equal(dim(pred4), c(6,3))
    expect_equal(colnames(pred4), c('mean', 's2', 'se'))

    # Check kernel$k matches when giving in as matrix or vector
    kn1 <- 5
    kn2 <- 7
    k_mat1 <- matrix(runif(kn1*d), nrow=kn1, ncol=d)
    k_mat2 <- matrix(runif(kn2*d), nrow=kn2, ncol=d)
    k_vec1 <- runif(d)
    k_vec2 <- runif(d)
    expect_equal(
      outer(1:kn1, 1:kn2, Vectorize(function(ii,jj) {gp$kernel$k(k_mat1[ii,], k_mat2[jj,])})),
      gp$kernel$k(k_mat1, k_mat2)
    )
    expect_equal(
      c(outer(1, 1:kn2, Vectorize(function(ii,jj) {gp$kernel$k(k_vec1, k_mat2[jj,])}))),
      c(gp$kernel$k(k_vec1, k_mat2))
    )
    expect_equal(
      c(outer(1:kn1, 1, Vectorize(function(ii,jj) {gp$kernel$k(k_mat1[ii,], k_vec2)}))),
      gp$kernel$k(k_mat1, k_vec2)
    )
    expect_equal(
      t(gp$kernel$k(k_mat2, k_mat1)),
      gp$kernel$k(k_mat1, k_mat2)
    )
    expect_equal(
      gp$kernel$k(k_vec2, k_mat1),
      gp$kernel$k(k_mat1, k_vec2)
    )


    # Check basics, but not for all kernels
    if (j<2.5) {
      expect_error(gp$plotLOO(), NA)
      expect_error(gp$plotmarginal(), NA)
      expect_error(gp$plotmarginalrandom(), NA)
      expect_error(gp$plot2D(), NA)
      expect_error(plot(gp), NA)
    }

    # Check EI for some kernels
    if (j<2.5) {
      expect_error(mei1 <- gp$maxEI(), NA)
      expect_is(mei1, "list")
      expect_equal(length(mei1), 2)
      expect_equal(length(mei1$par), 2)
      expect_equal(length(mei1$val), 1)
      expect_error(gp$maxqEI(npoints=1), NA)
      expect_error(mei2 <- gp$maxqEI(npoints=2), NA)
      expect_is(mei2, "list")
      expect_equal(length(mei2), 2)
      expect_is(mei2$par, "matrix")
      expect_equal(dim(mei2$par), c(2,2))
      expect_equal(length(mei2$val), 1)
    }

    # Check kernel
    expect_error({kernprint <- capture_output(print(gp$kernel))}, NA)
    expect_is(kernprint, 'character')

    df <- gp$deviance()
    dg <- gp$deviance_grad(nug.update = T)
    dfg <- gp$deviance_fngr(nug.update = T)
    expect_equal(df, dfg$fn)
    expect_equal(dg, dfg$gr, tolerance = 1e-4)
    # Now check numeric gradient
    # Nugget gradient
    eps <- 1e-6 # 1e-8 was too small, caused errors
    kernpars <- gp$kernel$param_optim_start(jitter=T)
    nuglog <- -3.3
    gpd1 <- gp$deviance(nuglog = nuglog - eps/2, params = kernpars)
    gpd2 <- gp$deviance(nuglog = nuglog + eps/2, params= kernpars)
    numder <- (gpd2-gpd1)/eps
    actder <- gp$deviance_grad(nuglog=nuglog, params=kernpars, nug.update = T)
    expect_equal(numder, actder[length(actder)], 1e-3)
    # kernel params
    # Use a random value since it's not exact enough to work at minimum
    kernpars <- gp$kernel$param_optim_start(jitter=T)
    if (kern_char=="RatQuad") {
      # Make sure kernpars are reasonable to avoid error: logalpha not too big
      # cat('ratquad kernpars are', kernpars, "\n")
    }
    actgrad <- gp$deviance_grad(params = kernpars, nug.update = F)
    for (i in 1:length(kernpars)) {
      epsvec <- rep(0, length(kernpars))
      epsvec[i] <- eps
      dp1 <- gp$deviance(params = kernpars - epsvec/2)
      dp2 <- gp$deviance(params = kernpars + epsvec/2)
      # numgrad <- (dp2-dp1) / eps
      # Switching to 4-point to reduce numerical errors, helps a lot
      dp3 <- gp$deviance(params = kernpars - epsvec)
      dp4 <- gp$deviance(params = kernpars + epsvec)
      numgrad <- (-dp4 + 8*dp2 - 8*dp1 + dp3)/(12*eps/2)
      # cat(j, kern_char, i, numgrad, actgrad[i+1], abs((numgrad - actgrad[1+i])/numgrad), "\n")
      expect_equal(numgrad, actgrad[1+i], tolerance = 1e-2,
                   label=paste(j,kern_char,i, 'numgrad'))
    }
  }
})

# Check factor kernels ----
test_that("check factor kernels alone", {
  n <- 20
  d <- 1
  x <- matrix(runif(n*d), ncol=d)
  # second is factor dim
  x[, 1] <- sample(1:3, n, T)
  # f <- function(x) {abs(sin(x[1]^.8*6))^1.2 + log(1+x[2]) + x[1]*x[2]}
  f <- function(x) {x[1]^.7}
  y <- apply(x, 1, f) + rnorm(n,0,1e-4) #f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  kern_chars <- c('FactorKernel', 'OrderedFactorKernel',
                  'LatentFactorKernel', 'LatentFactorKernel')
  kern_list <- list(
    FactorKernel$new(D=1, nlevels=3, xindex=1),
    OrderedFactorKernel$new(D=1, nlevels=3, xindex=1),
    LatentFactorKernel$new(D=1, nlevels=3, xindex=1, latentdim = 1, s2_est = F),
    LatentFactorKernel$new(D=1, nlevels=3, xindex=1, latentdim = 2)
  )
  for (j in 1:length(kern_chars)) {
    kern_char <- kern_chars[j]
    if (exists('printkern') && printkern) cat(j, kern_char, "\n")
    # kern <- eval(parse(text=kern_char))
    # expect_is(kern, "R6ClassGenerator")
    # kern1 <- IgnoreIndsKernel$new(Gaussian$new(D=1), ignoreinds = 2)
    kern <- kern_list[[j]]

    gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=kern, parallel=FALSE,
                                  verbose=0, nug.est=T, restarts=0)
    expect_is(gp, "GauPro")
    expect_is(gp, "R6")

    # Check kernel
    expect_error({kernprint <- capture_output(print(gp$kernel))}, NA)
    expect_is(kernprint, 'character')

    df <- gp$deviance()
    dg <- gp$deviance_grad(nug.update = T)
    dfg <- gp$deviance_fngr(nug.update = T)
    expect_equal(df, dfg$fn)
    expect_equal(dg, dfg$gr, tolerance = 1e-4)
    # Now check numeric gradient
    # Nugget gradient
    eps <- 1e-4 # 1e-8 was too small, caused errors
    kernpars <- gp$kernel$param_optim_start(jitter=T)
    nuglog <- -3.3
    gpd1 <- gp$deviance(nuglog = nuglog - eps/2, params = kernpars)
    gpd2 <- gp$deviance(nuglog = nuglog + eps/2, params= kernpars)
    numder <- (gpd2-gpd1)/eps
    actder <- gp$deviance_grad(nuglog=nuglog, params=kernpars, nug.update = T)
    expect_equal(numder, actder[length(actder)], 1e-3)
    # kernel params
    # Use a random value since it's not exact enough to work at minimum
    kernpars <- gp$kernel$param_optim_start(jitter=T)
    actgrad <- gp$deviance_grad(params = kernpars, nug.update = F)
    for (i in 1:length(kernpars)) {
      epsvec <- rep(0, length(kernpars))
      epsvec[i] <- eps
      dp1 <- gp$deviance(params = kernpars - epsvec/2)
      dp2 <- gp$deviance(params = kernpars + epsvec/2)
      # numgrad <- (dp2-dp1) / eps
      # Switching to 4-point to reduce numerical errors, helps a lot
      dp3 <- gp$deviance(params = kernpars - epsvec)
      dp4 <- gp$deviance(params = kernpars + epsvec)
      numgrad <- (-dp4 + 8*dp2 - 8*dp1 + dp3)/(12*eps/2)
      # cat(j, kern_char, i, numgrad, actgrad[i+1],
      #     abs((numgrad - actgrad[1+i])/numgrad), "\n")
      expect_equal(numgrad, actgrad[1+i], tolerance = 1e-2,
                   label=paste(j,kern_char,i, 'numgrad'))
      # debugonce(gp$kernel$dC_dparams)
    }
  }
})

# Check factor kernels in products ----
test_that("check factor kernels in product", {
  n <- 20
  d <- 2
  x <- matrix(runif(n*d), ncol=d)
  # second is factor dim
  nlev <- 2
  x[, 2] <- sample(1:nlev, n, T)
  f <- function(x) {abs(sin(x[1]^.8*6))^1.2 + log(1+x[2]) + x[1]*x[2]}
  y <- apply(x, 1, f) + rnorm(n,0,1e-4) #f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  kern_chars <- c('FactorKernel', 'OrderedFactorKernel',
                  'LatentFactorKernel1', 'LatentFactorKernel2')
  kern_list <- list(
    FactorKernel$new(D=2, nlevels=nlev, xindex=2),
    OrderedFactorKernel$new(D=2, nlevels=nlev, xindex=2),
    LatentFactorKernel$new(D=2, nlevels=nlev, xindex=2, latentdim = 1),
    LatentFactorKernel$new(D=2, nlevels=nlev, xindex=2, latentdim = 2)
  )
  for (j in 1:length(kern_chars)) {
    kern_char <- kern_chars[j]
    if (exists('printkern') && printkern) cat(j, kern_char, "\n")
    # kern <- eval(parse(text=kern_char))
    # expect_is(kern, "R6ClassGenerator")
    kern1 <- IgnoreIndsKernel$new(Matern32$new(D=1), ignoreinds = 2)
    kern <- kern1 * kern_list[[j]]
    # kern <- kern_list[[j]]

    gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=kern, parallel=FALSE,
                                  verbose=0, nug.est=T, restarts=0)
    expect_is(gp, "GauPro")
    expect_is(gp, "R6")

    # Check kernel
    expect_error({kernprint <- capture_output(print(gp$kernel))}, NA)
    expect_is(kernprint, 'character')

    # Check EI
    expect_error(mei1 <- gp$maxEI(), NA)
    # qEI just adds more time
    # expect_error(gp$maxqEI(npoints=1), NA)
    # expect_error(gp$maxqEI(npoints=2), NA)

    df <- gp$deviance()
    dg <- gp$deviance_grad(nug.update = T)
    dfg <- gp$deviance_fngr(nug.update = T)
    expect_equal(df, dfg$fn)
    expect_equal(dg, dfg$gr, tolerance = 1e-4)
    # Now check numeric gradient
    # Nugget gradient
    eps <- 1e-4 # 1e-8 was too small, caused errors
    kernpars <- gp$kernel$param_optim_start(jitter=T)
    nuglog <- -3.3
    gpd1 <- gp$deviance(nuglog = nuglog - eps/2, params = kernpars)
    gpd2 <- gp$deviance(nuglog = nuglog + eps/2, params= kernpars)
    numder <- (gpd2-gpd1)/eps
    actder <- gp$deviance_grad(nuglog=nuglog, params=kernpars, nug.update = T)
    expect_equal(numder, actder[length(actder)], 1e-3)
    # kernel params
    # Use a random value since it's not exact enough to work at minimum
    kernpars <- gp$kernel$param_optim_start(jitter=T)
    actgrad <- gp$deviance_grad(params = kernpars, nug.update = F)
    numgrads <- actgrad*NA
    for (i in 1:length(kernpars)) {
      epsvec <- rep(0, length(kernpars))
      epsvec[i] <- eps
      dp1 <- gp$deviance(params = kernpars - epsvec/2)
      dp2 <- gp$deviance(params = kernpars + epsvec/2)
      # numgrad <- (dp2-dp1) / eps
      # Switching to 4-point to reduce numerical errors, helps a lot
      dp3 <- gp$deviance(params = kernpars - epsvec)
      dp4 <- gp$deviance(params = kernpars + epsvec)
      numgrad <- (-dp4 + 8*dp2 - 8*dp1 + dp3)/(12*eps/2)
      numgrads[1+i] <- numgrad
      # cat(j, kern_char, i, numgrad, actgrad[i+1],
      #     abs((numgrad - actgrad[1+i])/numgrad), "\n")
      expect_equal(numgrad, actgrad[1+i], tolerance = 1e-2,
                   label=paste(j,kern_char,i, 'numgrad'))
    }
  }
})

# Check product kernels behave properly ----
test_that("check product kernels behave properly", {
  gk1 <- Gaussian$new(D=1)
  gk2 <- Gaussian$new(D=1)
  expect_true(gk1$s2_est)
  expect_true(gk2$s2_est)
  gk12 <- gk1*gk2
  expect_true(gk1$s2_est)
  # gk2 gets flipped
  expect_true(!gk2$s2_est)
  expect_true(gk12$s2_est)

  # Now turn off s2_est on gk12
  gk12$s2_est <- F
  expect_true(!gk1$s2_est)
  expect_true(!gk2$s2_est)
  expect_true(!gk12$s2_est)

  # Turn back on
  gk12$s2_est <- TRUE
  expect_true(gk1$s2_est)
  expect_true(!gk2$s2_est)
  expect_true(gk12$s2_est)


  h1 <- IgnoreIndsKernel$new(Gaussian$new(D=2), 1)
  h2 <- Matern32$new(D=1)
  h3 <- PowerExp$new(D=1)
  expect_true(h1$s2_est)
  expect_true(h2$s2_est)
  expect_true(h3$s2_est)
  h <- h1*h2*h3
  expect_true(h1$s2_est)
  expect_true(!h2$s2_est)
  expect_true(!h3$s2_est)
  expect_true(h$s2_est)

  # Turn off
  h$s2_est <- F
  expect_true(!h1$s2_est)
  expect_true(!h2$s2_est)
  expect_true(!h3$s2_est)
  expect_true(!h$s2_est)

  # Turn on
  h$s2_est <- T
  expect_true(h1$s2_est)
  expect_true(!h2$s2_est)
  expect_true(!h3$s2_est)
  expect_true(h$s2_est)
})

# Formula/data input ----
test_that("Formula/data input", {
  n <- 30
  tdf <- data.frame(a=runif(n), b=runif(n), c=factor(sample(5:6,n,T)), d=runif(n), e=sample(letters[1:3], n,T))
  tdf$z <- with(tdf, a+a*b+b^2)
  gpf <- GauPro_kernel_model$new(X=tdf, Z=z ~ a + b + c + e, kernel='gauss')
  expect_true("GauPro" %in% class(gpf))
  expect_equal(ncol(gpf$X), 4)
  expect_true(is.matrix(gpf$X))
  expect_error(predict(gpf, tdf), NA)
})

# EI ----
test_that("EI with mixopt", {
  n <- 30
  tdf <- data.frame(a=runif(n), b=runif(n, -1,1),
                    c=(sample(5:6,n,T)),
                    d=sample(c(.1,.2,.3,.4), n, T),
                    # e=sample(letters[1:3], n,T),
                    f=sample(10:30, n, T))
  z <- with(tdf, a+a*b+b^2 +5*(d-.22)^2*(f-22)^2)
  gpf <- GauPro_kernel_model$new(X=as.matrix(tdf), Z=z, kernel='m52')
  gpf$maxEI(minimize = T)
  mop <- c(
    mixopt::mopar_cts(0,1),
    mixopt::mopar_cts(-1,1),
    mixopt::mopar_unordered(5:6),
    mixopt::mopar_ordered(c(.1,.2,.3,.4)),
    mixopt::mopar_ordered(10:30)
  )
  expect_error(gpf$maxEI(mopar = mop, minimize = T), NA)
  expect_error(gpf$maxqEI(npoints = 3, mopar = mop, minimize = T), NA)
})
