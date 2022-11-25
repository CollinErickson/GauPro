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
                    Gaussian$new(D=2) * PowerExp$new(D=2),
                    Matern52$new(D=2) + Matern32$new(D=2))
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

    expect_no_warning({
      expect_error({
        gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=kern, parallel=FALSE,
                                      verbose=0, nug.est=T, restarts=0)
      }, NA)
    })
    expect_is(gp, "GauPro")
    expect_is(gp, "R6")

    # Check predict
    expect_error(pred1 <- predict(gp, runif(2)), NA)
    expect_true(is.numeric(pred1))
    expect_true(!is.matrix(pred1))
    expect_equal(length(pred1), 1)
    # Predict with SE
    expect_error(pred2 <- predict(gp, runif(2), se.fit=T), NA)
    expect_true(is.data.frame(pred2))
    expect_equal(dim(pred2), c(1,3))
    expect_equal(colnames(pred2), c('mean', 's2', 'se'))
    # Matrix
    expect_error(pred3 <- predict(gp, matrix(runif(12), ncol=2)), NA)
    expect_true(is.numeric(pred3))
    expect_true(!is.matrix(pred3))
    expect_equal(length(pred3), 6)
    # Matrix with SE
    expect_error(pred4 <- predict(gp, matrix(runif(12), ncol=2), se.fit=T), NA)
    expect_true(is.data.frame(pred4))
    expect_equal(dim(pred4), c(6,3))
    expect_equal(colnames(pred4), c('mean', 's2', 'se'))
    # Predict mean dist
    expect_error(pred5 <- predict(gp, matrix(runif(12), ncol=2), se.fit=T,
                                  mean_dist=TRUE), NA)

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

    # Summary
    expect_no_warning(
      expect_error(capture.output(summary(gp)), NA)
    )

    # Kernel plot
    expect_error(plot(gp$kernel), NA)

    # Test importance
    expect_error(capture.output(imp <- gp$importance(plot=F)), NA)
    expect_true(is.numeric(imp))
    expect_equal(names(imp), c("X1", "X2"))

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

    # Test grad. Implicitly tests kernel$dC_dx.
    if (j %in% c(1:3, 6:7, 9, 11:13)) {
      xgrad <- runif(2) #matrix(runif(6), ncol=2)
      expect_no_error(symgrad <- gp$grad(xgrad))
      expect_equal(numDeriv::grad(gp$pred, x=xgrad),
                   c(symgrad),
                   tolerance=1e-4)
      # grad at self shouldn't be zero, except for Exponential
      expect_no_error(gpgradX <- gp$grad(gp$X))
      if (j != 9) {
        expect_true(!any(is.na(gpgradX)))
      }
    } else {
      if (exists('printkern') && printkern) {
        cat("grad/dC_dx not tested for", j, kern_char, "\n")
      }
    }

    # Test gradpredvar
    # FIX m32/m52
    if (j %in% c(1:1, 6:7, 9, 11:12)) {
      expect_no_error(gpv <- gp$gradpredvar(xgrad))
      # numDeriv::grad(func=function(x) gp$pred(x, se=T)$s2, gpv)
      npv <- 39
      XXpv <- matrix(runif(2*npv), ncol=2)
      expect_no_error(XXgpv <- gp$gradpredvar(XXpv))
      gpvmatches <- 0
      for (iii in 1:npv) {
        numpv <- numDeriv::grad(func=function(x) {gp$pred(x, se=T)$s2},
                                XXpv[iii,])
        # expect_equal(numpv, XXgpv[iii,], tolerance = 1e-2)
        pvclose <- all(ifelse(XXgpv[iii,]==0,
                              abs(numpv - XXgpv[iii,]) < 1e-8,
                              abs((numpv - XXgpv[iii,]) / XXgpv[iii,]) < 1e-2))
        if (pvclose) {gpvmatches <- gpvmatches + 1}
      }
      if (exists('printkern') && printkern) {
        cat(kern_char, "gpv close on ", gpvmatches, "/", npv, "\n")
      }
      expect_true(gpvmatches > npv/2)
    } else {
      if (exists('printkern') && printkern) {
        cat("gradpredvar not tested for", j, kern_char, "\n")
      }
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

    # max attempts
    maxattempts <- 10
    numgradtol <- 1e-3
    for (iatt in 1:maxattempts) {
      goodsofar <- TRUE
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
        # expect_equal(numgrad, actgrad[1+i], tolerance = 1e-6,
        # label=paste(j,kern_char,i, 'numgrad'))
        alleq <- all.equal(actgrad[1+i], numgrad, tolerance=numgradtol)
        if (isTRUE(alleq)) {
          expect_equal(numgrad, actgrad[1+i], tolerance = numgradtol,
                       label=paste(j,kern_char,i, 'numgrad'))
        } else {
          if (exists('printkern') && printkern) {
            cat("FAILURE", kern_char, iatt, i, alleq, "\n")
          }
          goodsofar <- FALSE
        }
      }
      if (goodsofar) {
        break
      } else {
        if (exists("printkern") && printkern) {
          cat("failed on", kern_char, "attempt", iatt,
              alleq,
              (numgrad-actgrad[1+i]) / actgrad[1+i],
              "\n")
        }
      }
    }
    # This is where it will fail if none succeeded
    expect_true(alleq,
                label=paste(kern_char,
                            'numgrad matches symbolic grad (failed on all',
                            maxattempts, "attempts)"))
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
  y <- apply(x, 1, f) + rnorm(n,0,1e-1) #f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  kern_chars <- c('FactorKernel', 'OrderedFactorKernel',
                  'LatentFactorKernel', 'LatentFactorKernel')
  kern_list <- list(
    FactorKernel$new(D=1, nlevels=3, xindex=1),
    OrderedFactorKernel$new(D=1, nlevels=3, xindex=1),
    LatentFactorKernel$new(D=1, nlevels=3, xindex=1, latentdim = 1, s2_est=F, s2=.3),
    LatentFactorKernel$new(D=1, nlevels=3, xindex=1, latentdim = 2)
  )
  for (j in 1:length(kern_chars)) {
    kern_char <- kern_chars[j]
    if (exists('printkern') && printkern) cat(j, kern_char, "\n")
    # kern <- eval(parse(text=kern_char))
    # expect_is(kern, "R6ClassGenerator")
    # kern1 <- IgnoreIndsKernel$new(Gaussian$new(D=1), ignoreinds = 2)
    kern <- kern_list[[j]]

    expect_error({
      gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=kern, parallel=FALSE,
                                    verbose=0, nug.est=T, restarts=0,
                                    track_optim = !T)
    }, NA)
    # gp$plot1D()
    # gp$plot_track_optim()
    expect_is(gp, "GauPro")
    expect_is(gp, "R6")

    # Basic check for k
    expect_equal(kern$k(1, 1), c(kern$k(matrix(1), matrix(1))), tolerance=1e-4)

    # Test predict
    if (T || (j %in% 1:4)) {
      expect_error(predict(gp, 1:3, se.fit = T), NA,
                   info = paste("bad pred in", kern_char))
      expect_warning(predict(gp, 1:3, se.fit = T), NA,
                     info = paste("bad pred in", kern_char))
      # Test plot
      expect_no_error(pp <- gp$plot1D())
      expect_no_error(suppressMessages({pp <- plot(gp)}))
    } else {
      if (exists('printkern') && printkern) {
        cat("factorkernel", j, kern_char, "not testing pred/plot", "\n")
      }
    }

    # Test importance
    expect_error(capture.output(imp <- gp$importance(plot=F)), NA)
    expect_true(is.numeric(imp))
    expect_equal(names(imp), c("X1"))

    # Summary works and prints
    expect_no_warning(
      expect_error(capture.output(summary(gp)), NA)
    )

    # grad should be NA
    expect_no_error(grd <- gp$grad(gp$X))
    expect_true(is.matrix(grd))
    expect_equal(dim(grd), c(n, 1))
    expect_true(all(is.na(grd)))

    # Check kernel
    expect_error({kernprint <- capture_output(print(gp$kernel))}, NA)
    expect_is(kernprint, 'character')
    expect_no_error(plot(gp$kernel))

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

    expect_no_error({
      gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=kern, parallel=FALSE,
                                    verbose=0, nug.est=T, restarts=0)
    })
    expect_is(gp, "GauPro")
    expect_is(gp, "R6")

    # Check kernel print
    expect_error({kernprint <- capture_output(print(gp$kernel))}, NA)
    expect_is(kernprint, 'character')

    # Summary
    expect_no_warning(
      expect_error(capture.output(summary(gp)), NA)
    )

    # Check LOO
    expect_no_error(capture.output(gp$plotLOO()))

    # Check EI
    expect_error(mei1 <- gp$maxEI(), NA)
    # qEI just adds more time
    # expect_error(gp$maxqEI(npoints=1), NA)
    # expect_error(gp$maxqEI(npoints=2), NA)

    # grad should be NA on factor column, but not other
    expect_no_error(grd <- gp$grad(cbind(runif(10), sample(1:2, 10, T))))
    expect_true(is.matrix(grd))
    expect_equal(dim(grd), c(10, d))
    expect_true(all(is.na(grd[, 2])))
    expect_true(all(!is.na(grd[, 1])))

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
  tdf <- data.frame(a=runif(n), b=runif(n), c=factor(sample(5:6,n,T)),
                    d=runif(n), e=sample(letters[1:3], n,T))
  tdf$z <- with(tdf, a+a*b+b^2)
  gpf <- GauPro_kernel_model$new(X=tdf, Z=z ~ a + b + c + e, kernel='gauss')
  expect_true("GauPro" %in% class(gpf))
  expect_equal(ncol(gpf$X), 4)
  expect_true(is.matrix(gpf$X))
  expect_error(predict(gpf, tdf), NA)
})

test_that("Formula/data input 2", {
  library(dplyr)
  n <- 63
  xdf <- tibble(
    a=rnorm(n),
    b=runif(n),
    c=sample(letters[1:5], n, T),
    d=factor(sample(letters[6:9], n, T)),
    e=rexp(n),
    # f=rnorm(n),
    z=a*b + a^2*ifelse(c %in% c('a', 'b'), 1, .5) +
      e*ifelse(d %in% c('g','h'), 1, -1) +
      ifelse(paste0(d,e) %in% c('af', 'ah', 'cf', 'cg', 'ci'),4,0) +
      rnorm(n, 1e-3)
  )
  xdf
  # xdf %>% str
  # xdf %>% GGally::ggpairs()

  # Test fit
  expect_error(gpdf <- GauPro_kernel_model$new(z ~ ., data=xdf, kernel='m32'), NA)
  expect_true("formula" %in% class(gpdf$formula))
  expect_true("terms" %in% class(gpdf$formula))
  # Kernel should automatically have factors for chars
  expect_true("GauPro_kernel_product" %in% class(gpdf$kernel), NA)
  # Test predict
  expect_error(predict(gpdf, xdf), NA)
  # Missing col, give descriptive error
  expect_error(predict(gpdf, xdf[,1:3]))
  # Test pred LOO
  expect_error(gpdf$plotLOO(), NA)
  # Test importance
  expect_error(capture.output(imp <- gpdf$importance(plot=F)), NA)
  expect_true(is.numeric(imp))
  expect_equal(names(imp), attr(gpdf$formula, "term.labels"))
  # Summary
  expect_no_warning(
    expect_error(capture.output(summary(gpdf)), NA)
  )
  # Test EI
  expect_no_error(gpdf$EI(xdf[1,]))
  # Test maxEI
  expect_error(dfEI <- gpdf$maxEI(), NA)
  expect_true(is.data.frame(dfEI$par))
  expect_equal(colnames(dfEI$par), colnames(xdf)[1:5])
  expect_equal(dim(dfEI$par), c(1,5))
  rm(dfEI)
  # maxEI with mopar
  expect_error(dfEI <- gpdf$maxEI(
    mopar = c(mixopt::mopar_cts(-3,3),
              mixopt::mopar_cts(0,1),
              mixopt::mopar_unordered(letters[1:5]),
              mixopt::mopar_unordered(letters[6:9]),
              mixopt::mopar_cts(0,4)
    )
  ), NA)
  expect_true(is.data.frame(dfEI$par))
  expect_equal(colnames(dfEI$par), colnames(xdf)[1:5])
  expect_equal(dim(dfEI$par), c(1,5))
  # Test qEI with mopar
  # expect_error(dfqEI <- gpdf$maxqEI(
  #   npoints = 2,
  #   mopar = c(mixopt::mopar_cts(-3,3),
  #             mixopt::mopar_cts(0,1),
  #             mixopt::mopar_unordered(letters[1:5]),
  #             mixopt::mopar_unordered(letters[6:9]),
  #             mixopt::mopar_cts(0,4)
  #   )
  # ), NA)
  # expect_true(is.data.frame(dfqEI$par))
  # expect_equal(colnames(dfqEI$par), colnames(xdf)[1:5])
  # expect_equal(dim(dfqEI$par), c(2,5))


  # Try other arg names
  expect_error(gpdf <- GauPro_kernel_model$new(z ~ ., xdf, kernel='m32'), NA)
  expect_error(gpdf <- GauPro_kernel_model$new(xdf, z ~ ., kernel='m32'), NA)
  expect_error(gpdf <- GauPro_kernel_model$new(formula=z ~ ., data=xdf, kernel='m32'), NA)
  expect_error(gpdf <- GauPro_kernel_model$new(z ~ ., xdf, kernel='m32'), NA)
  expect_error(gpdf <- GauPro_kernel_model$new(z ~ xdf$a + xdf$c + xdf$e, xdf, kernel='m32'), NA)
  expect_error(gpdf <- GauPro_kernel_model$new(z ~ ., data=xdf, kernel='m32'), NA)
  rm(gpdf, dfEI)

  # Only fit on chars
  expect_error(gpch <- GauPro_kernel_model$new(z ~ c + d, data=xdf, kernel='m32'), NA)
  expect_error(gpch <- GauPro_kernel_model$new(z ~ a + c, data=xdf, kernel='m32'), NA)

  # Try more complex
  # Give error if formula includes
  expect_error(gpco <- GauPro_kernel_model$new(z ~ a*b + c + d, data=xdf, kernel='m32'))
  expect_error(gpco <- GauPro_kernel_model$new(z ~ exp(a) + b*e + c, data=xdf, kernel='m32'))
  # Transformations work fine
  expect_error(gpco <- GauPro_kernel_model$new(exp(z) ~ exp(a) + c, data=xdf, kernel='m32'), NA)
  expect_equal(exp(xdf$z), c(gpco$Z))
  expect_equal(exp(xdf$a), unname(gpco$X[,1]))
  expect_error(gpco <- GauPro_kernel_model$new(exp(z) ~ exp(a) + c +sqrt(e), data=xdf, kernel='m32'), NA)

  # Transf
  expect_error(gpdf <- GauPro_kernel_model$new(abs(z) ~ exp(a) + c + d + sqrt(e), data=xdf, kernel='m32'), NA)
  # Test EI
  expect_error(dfEI <- gpdf$maxEI(), NA)
  expect_true(is.data.frame(dfEI$par))
  expect_equal(colnames(dfEI$par), attr(gpdf$formula, "term.labels"))
  expect_equal(dim(dfEI$par), c(1,4))
  # Test qEI
  # expect_error(dfqEI <- gpdf$maxqEI(npoints = 2), NA)
  # expect_true(is.data.frame(dfqEI$par))
  # expect_equal(colnames(dfqEI$par), attr(gpdf$formula, "term.labels"))
  # expect_equal(dim(dfqEI$par), c(2,4))
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

# Misc ----
test_that("S3 +/*", {
  w1 <- Cubic$new(D=4)
  expect_error(w1 * 1, NA)
  expect_error(w1 * 2)
  expect_error(w1 + 0, NA)
  expect_error(w1 + 1)
})
test_that("IgnoreInds", {
  r1 <- Cubic$new(D=1)
  # Ignore inds must be vector of numeric
  expect_error(r2 <- IgnoreIndsKernel$new(k=r1, ignoreinds = list(2)))
  expect_error(IgnoreIndsKernel$new(k=r1, ignoreinds=1.00001))
})
test_that("Wide range", {
  d <- 3
  n <- 40
  x <- cbind(runif(n, 50, 150),
             runif(n, -9,-7),
             runif(n, 1200, 3000))
  f <- function(x) {x[1] + 14* x[2]^2 + x[3]}
  y <- apply(x, 1, f) + rnorm(n, 0, 1)
  expect_no_error(e1 <- GauPro_kernel_model$new(
    x, y, kernel=Gaussian$new(D=3, s2=3e7, s2_lower=1e7),
    verbose=0, restarts=25))
  # e1$plotLOO(); print(e1$s2_hat); print(e1$nug)
  expect_no_error(e1$summary())
  expect_no_error(e1$plotLOO())
  expect_no_error(e1$plotmarginalrandom())
  expect_true(mean(abs((e1$Z-e1$pred(e1$X)) / e1$Z)) < 1e-2)
})
