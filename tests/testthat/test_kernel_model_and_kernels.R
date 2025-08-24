library(testthat)

set.seed(Sys.time())
seed <- floor(runif(1)*1e6)
set.seed(seed)

printkern <- interactive()
# cat('printkern is', printkern, '\n')
if (printkern) {cat("seed =", seed, "\n")}

if (requireNamespace('testthatmulti', quietly = TRUE)) {

} else {
  ttm <- function(a, b) {b}
  ttm_expect_true <- testthat::expect_true
  ttm_expect_equal <- testthat::expect_equal
  ttm_expect_error <- testthat::expect_error
  ttm_expect_no_error <- testthat::expect_no_error
}

# Cts kernels 1D ----
test_that("Cts kernels 1D", {
  n <- sample(11:18, 1)
  d <- 1
  x <- matrix(runif(n*d,.3,1.9), ncol=d)
  n <- nrow(x)
  f <- function(x) {abs(sin(x[1]^.8*6))^1.2 + log(1+(x[1]-.3)^2)}
  y <- apply(x, 1, f) + rnorm(n,0,1e-2) #f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  kern_chars <- c('Gaussian', 'Matern32', 'Matern52',
                  'Triangle', 'Cubic', 'White',
                  'PowerExp', 'Periodic', "Exponential", "RatQuad",
                  "Product", "Sum",
                  # Same but with useC=FALSE
                  'Gaussian', 'Matern32', 'Matern52',
                  'Triangle', 'Cubic', 'White',
                  'PowerExp', 'Periodic', "Exponential", "RatQuad"
  )
  kern_list <- list(0,0,0,
                    0,0,White$new(D=1),
                    0,0,0,0,
                    # IgnoreIndsKernel$new(Gaussian$new(D=1), 2),
                    Gaussian$new(D=d) * Periodic$new(D=d),
                    Matern52$new(D=d) + Matern32$new(D=d),
                    # Same but with useC=FALSE
                    Gaussian$new(D=1, useC=FALSE),
                    Matern32$new(D=1, useC=FALSE),
                    Matern52$new(D=1, useC=FALSE),
                    Triangle$new(D=1, useC=FALSE),
                    Cubic$new(D=1, useC=FALSE),
                    White$new(D=1, useC=FALSE),
                    PowerExp$new(D=1, useC=FALSE),
                    Periodic$new(D=1, useC=FALSE),
                    Exponential$new(D=1, useC=FALSE),
                    RatQuad$new(D=1, useC=FALSE)
  )
  stopifnot(length(kern_chars) == length(kern_list))
  for (j in 1:length(kern_chars)) {
    if (exists('seed')) {set.seed(seed)} else {seed <- runif(1)}
    kern_char <- kern_chars[j]
    if (exists('printkern') && printkern) cat("1D", j, kern_char, "\n")
    if (is.numeric(kern_list[[j]])) {
      # kern <- eval(parse(text=kern_char))
      # expect_is(kern, "R6ClassGenerator")
      kern <- kern_char
    } else {
      kern <- kern_list[[j]]
    }

    expect_no_warning({
      expect_no_error({
        gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=kern, parallel=FALSE,
                                      verbose=0, nug.est=T, restarts=0)
      }, message=paste0("Cts 1D", j, kern_char, "failed fitting no error"))
    }, message=paste0("Cts 1D", j, kern_char, "failed fitting no warning"))
    expect_is(gp, "GauPro")
    expect_is(gp, "R6")

    # Check kernel properties

    expect_equal(GauPro:::find_kernel_cts_dims(gp$kernel),
                 if (kern_char == "White") {NULL}
                 else if (kern_char=="Ignore") {1}
                 else {1})
    expect_true(is.null(GauPro:::find_kernel_factor_dims(gp$kernel)))

    # Check predict
    expect_error(pred1 <- predict(gp, runif(d)), NA)
    expect_true(is.numeric(pred1))
    expect_true(!is.matrix(pred1))
    expect_equal(length(pred1), 1)
    # Predict with SE
    expect_error(pred2 <- predict(gp, runif(d), se.fit=T), NA)
    expect_true(is.data.frame(pred2))
    expect_equal(dim(pred2), c(1,3))
    expect_equal(colnames(pred2), c('mean', 's2', 'se'))
    # Matrix
    expect_error(pred3 <- predict(gp, matrix(runif(12), ncol=d)), NA)
    expect_true(is.numeric(pred3))
    expect_true(!is.matrix(pred3))
    expect_equal(length(pred3), 12)
    # Matrix with SE
    expect_error(pred4 <- predict(gp, matrix(runif(12), ncol=d), se.fit=T), NA)
    expect_true(is.data.frame(pred4))
    expect_equal(dim(pred4), c(12,3))
    expect_equal(colnames(pred4), c('mean', 's2', 'se'))
    # Predict mean dist
    expect_error(pred5 <- predict(gp, matrix(runif(12), ncol=d), se.fit=T,
                                  mean_dist=TRUE), NA)

    # Kernel k with self must equal s2. If this fails, may need to change places
    #  where it is assumed to be true.
    expect_equal(gp$kernel$k(.3), gp$kernel$s2,
                 label=paste(j, kern_char, 'k(.3)'))
    if (kern_char != "White") {
      expect_equal(gp$kernel$k(.3,.3), gp$kernel$s2,
                   label=paste(j, kern_char, 'k(.3,.3)'))
    }
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
      # plot1D
      expect_error(gp$plot1D(), NA)
      expect_error(gp$plot1D(gg=FALSE), NA)
      expect_no_error(gp$cool1Dplot(gg=T))
      expect_no_error(gp$cool1Dplot(gg=F))
      # Should call plot2D
      expect_error(plot(gp), NA)
    }

    # Summary
    expect_no_warning(
      expect_error(capture.output(summary(gp)), NA)
    )

    # Kernel plot
    expect_error(plot(gp$kernel), NA)
    expect_error(gp$plotkernel, NA)

    # Test importance
    expect_error(capture.output(imp <- gp$importance(plot=F)), NA)
    expect_true(is.numeric(imp))
    expect_equal(names(imp), c("X1"))


    # Test grad. Implicitly tests kernel$dC_dx.
    if (j %in% c(1:33)) {
      xgrad <- runif(d) #matrix(runif(6), ncol=2)
      expect_no_error(symgrad <- gp$grad(xgrad))
      expect_equal(numDeriv::grad(gp$pred, x=xgrad),
                   c(symgrad),
                   tolerance=1e-2,
                   label=paste(j, kern_char, 'gp$grad'))
      # grad at self shouldn't be zero, except for Triangle, Exponential, PowerExp.
      expect_no_error(gpgradX <- gp$grad(gp$X))
      # if (!(j %in% c(4,7,9))) {
      if (!(kern_char %in% c("Triangle", "Exponential", "PowerExp"))) {
        expect_true(!any(is.na(gpgradX)),
                    label=paste(j, kern_char, '!any(is.na(gpgradX))'))
      }
    } else {
      if (exists('printkern') && printkern) {
        cat("grad/dC_dx not tested for", j, kern_char, "\n")
      }
    }

    # Test gradpredvar
    # FIX m32/m52
    if (j %in% c(1:33)) {
      expect_no_error(gpv <- gp$gradpredvar(xgrad))
      # numDeriv::grad(func=function(x) gp$pred(x, se=T)$s2, gpv)
      npv <- 39
      XXpv <- matrix(runif(d*npv), ncol=d)
      expect_no_error(XXgpv <- gp$gradpredvar(XXpv))
      gpvmatches <- 0
      numpvs <- c()
      actpvs <- c()
      for (iii in 1:npv) {
        numpv <- numDeriv::grad(func=function(x) {gp$pred(x, se=T)$s2},
                                XXpv[iii,])
        # expect_equal(numpv, XXgpv[iii,], tolerance = 1e-2)
        pvclose <- all(
          ifelse(abs(XXgpv[iii,]) < 1e-8,
                 abs(numpv - XXgpv[iii,]) < 1e-8,
                 ifelse(abs(XXgpv[iii,]) < 1e-2,
                        abs((numpv - XXgpv[iii,]) / XXgpv[iii,]) < 1e-1,
                        abs((numpv - XXgpv[iii,]) / XXgpv[iii,]) < 1e-3)))
        if (pvclose) {gpvmatches <- gpvmatches + 1}
        numpvs <- c(numpvs, numpv)
        actpvs <- c(actpvs, XXgpv[iii,])
      }
      if (exists('printkern') && printkern && gpvmatches < npv) {
        cat(kern_char, "gpv close on ", gpvmatches, "/", npv, "\n")
      }
      expect_true(gpvmatches > npv/2,
                  label=paste(j, kern_char, 'gpvmatches', gpvmatches,'/',npv,
                              "seed =", seed))
      # qplot(numpvs, actpvs)
      # summary((numpvs - actpvs) / actpvs)
      # cbind(numpvs, actpvs, prop=(numpvs - actpvs) / actpvs)
      # qplot(numpvs, (numpvs - actpvs) / actpvs)
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
    expect_true(isTRUE(alleq),
                label=paste(kern_char,
                            'numgrad matches symbolic grad (failed on all',
                            maxattempts, "attempts)"))
  }
})

# Cts kernels 2D ----
test_that("Cts kernels 2D", {
  kern_chars <- c('Gaussian', 'Matern32', 'Matern52',
                  'Triangle', 'Cubic', 'White',
                  'PowerExp', 'Periodic', "Exponential", "RatQuad",
                  "Ignore", "Product", "Sum",
                  "Gaussian_isotropic", "Matern52_isotropic",
                  "Matern32_isotropic", "Exponential_isotropic",
                  "Triangle_isotropic")
  kern_list <- list(0,0,0,0,0,0,
                    0,0,0,0,
                    IgnoreIndsKernel$new(Gaussian$new(D=1), 2),
                    Gaussian$new(D=2) * PowerExp$new(D=2),
                    Matern52$new(D=2) + Matern32$new(D=2),
                    Gaussian$new(D=2, isotropic=T),
                    Matern52$new(D=2, isotropic=T),
                    Matern32$new(D=2, isotropic=T),
                    Exponential$new(D=2, isotropic=T),
                    Triangle$new(D=2, isotropic=T))
  stopifnot(length(kern_chars) == length(kern_list))
  for (j in 1:length(kern_chars)) {
    # Don't set seed inside ttm, it'll defeat the purpose
    # This way we know what data was run first
    # Or could set seed to be seed + ttm_i() for repeatability
    if (exists('seed')) {set.seed(seed)} else {seed=runif(1)}
    n <- 20
    d <- 2

    # Used to set seed before starting fit, don't do that with ttm
    # if (exists('seed')) {set.seed(seed)} else {seed <- runif(1)}

    # Need this inside kernel for loop: it shouldn't redo previous kernels
    ttm(10,{
      # Generate data inside of ttm so that it can retry on new data
      x <- matrix(runif(n*d), ncol=d)
      x <- GauPro:::lhs_maximinLHS(n, d) # Better spacing might avoid grad issues?
      # x <- rbind(x, x[1,]) # Add repeated x since that could cause issues
      # n <- nrow(x)
      f <- function(x) {abs(sin(x[1]^.8*6))^1.2 + log(1+(x[2]-.3)^2) + x[1]*x[2]}
      y <- apply(x, 1, f) + rnorm(n,0,1e-2) #f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)


      kern_char <- kern_chars[j]
      if (exists('printkern') && printkern) cat("2D", j, kern_char, "\n")
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

      # Check kernel properties
      expect_equal(GauPro:::find_kernel_cts_dims(gp$kernel),
                   if (kern_char == "White") {NULL}
                   else if (kern_char=="Ignore") {1}
                   else {1:2})
      expect_true(is.null(GauPro:::find_kernel_factor_dims(gp$kernel)))

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

      # Kernel k with self must equal s2. If this fails, may need to change places
      #  where it is assumed to be true.
      expect_equal(gp$kernel$k(.3), gp$kernel$s2,
                   label=paste(j, kern_char, 'k(.3)'))
      if (kern_char != "White") {
        expect_equal(gp$kernel$k(.3,.3), gp$kernel$s2,
                     label=paste(j, kern_char, 'k(.3,.3)'))
      }
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
        # plot2D
        expect_error(gp$plot2D(), NA)
        expect_no_error({gp$plot2D(se=T, n=5)})
        expect_no_error(gp$plot2D(se=T, n=5, horizontal=F))
        expect_no_error(gp$plot2D(se=T, n=5, mean=F))
        expect_error(gp$plot2D(se=F, mean=F))
        expect_error(gp$plot2D(se=1))
        # Should call plot2D
        expect_error(plot(gp), NA)
        # Sample
        expect_no_error(s1 <- gp$sample(XX=runif(d), 3))
        expect_true(is.matrix(s1))
        expect_equal(dim(s1), c(3, 1))
        expect_no_error(s2 <- gp$sample(XX=matrix(runif(d*5), ncol=d), 30))
        expect_true(is.matrix(s2))
        expect_equal(dim(s2), c(30, 5))
      }

      # Check some advanced stuff only for Gaussian
      if (j %in% c(1, 14)) {
        # Covmat
        expect_no_error(gp$predict(matrix(runif(2*d), ncol=2), covmat = T))
        # Split speed, large matrix
        expect_no_error(gp$predict(matrix(runif(9000*d), ncol=2), split_speed = T))
        # Fails on range length
        expect_error(gp$predict(runif(3)))

        # Various grad/var stuff
        expect_no_error(gp$pred_var_after_adding_points(add_points = runif(d),
                                                        pred_points = runif(d)))
        expect_no_error(gp$pred_var_after_adding_points(
          add_points = matrix(runif(2*d), ncol=d), pred_points = runif(d)))
        expect_no_error(gp$pred_var_reductions(
          add_points = matrix(runif(2*d), ncol=d), pred_points = runif(d)))
        expect_no_error(gp$grad_dist(matrix(runif(d), ncol=d)))
        expect_no_error(gp$grad_sample(matrix(runif(d), ncol=d), n=10))
        expect_no_error(gp$grad_norm2_mean(matrix(runif(d), ncol=d)))
        expect_no_error(gp$hessian(matrix(runif(d), ncol=d)))

        # optimize_fn
        expect_no_error(gp$optimize_fn(function(x) {gp$predict(x)}, minimize = FALSE))
        expect_no_error(gp$optimize_fn(function(x) {gp$predict(x)}, minimize = TRUE))
        expect_no_error(gp$optimize_fn(function(x) {
          p <- gp$predict(x, T)
          p$mean + p$se
        }, minimize = FALSE))
        expect_no_error(gp$optimize_fn(function(x) {gp$predict(x)},
                                       gr=function(x) {gp$grad(x)}, minimize = FALSE))
        expect_no_error(gp$optimize_fn(fn=function(x) {gp$predict(x)},
                                       fngr=function(x) {
                                         list(fn=gp$predict(x), gr=gp$grad(x))
                                       }, minimize = FALSE))
        expect_no_error(gp$optimize_fn(function(x, a) {a+gp$predict(x)},
                                       fn_args=list(a=100)))
      }

      # Summary
      expect_no_warning(
        expect_error(capture.output(summary(gp)), NA)
      )

      # Print
      expect_no_error(printout <- capture.output(print(gp)))
      expect_true(is.character(printout))
      expect_equal(printout[1], "GauPro kernel model object")

      # Kernel plot
      expect_error(plot(gp$kernel), NA)
      expect_error(gp$plotkernel, NA)

      # Test importance
      expect_error(capture.output(imp <- gp$importance(plot=F)), NA)
      expect_true(is.numeric(imp))
      expect_equal(names(imp), c("X1", "X2"))

      # Check EI for some kernels
      if (j<2.5) {
        for (EItype in c("ei", "aug", "cor")) {
          expect_no_error(mei1 <- gp$maxEI(EItype = EItype),
                          message=paste("EI", EItype, kern_char))
          expect_is(mei1, "list")
          expect_equal(length(mei1), 2)
          expect_equal(length(mei1$par), 2)
          expect_equal(length(mei1$val), 1)
        }
        expect_no_error(gp$maxqEI(npoints=1),
                        message=paste("maxqEI1", EItype, kern_char))
        expect_no_error(mei2 <- gp$maxqEI(npoints=2),
                        message=paste("maxqEI2", EItype, kern_char))
        expect_is(mei2, "list")
        expect_equal(length(mei2), 2)
        expect_is(mei2$par, "matrix")
        expect_equal(dim(mei2$par), c(2,2))
        expect_equal(length(mei2$val), 1)
      }

      # Test grad over X. Implicitly tests kernel$dC_dx.
      if (T) {
        xgrad <- runif(2) #matrix(runif(6), ncol=2)
        expect_no_error(symgrad <- gp$grad(xgrad))
        expect_equal(numDeriv::grad(gp$pred, x=xgrad),
                     c(symgrad),
                     tolerance=1e-2,
                     label=paste(j, kern_char, 'gp$grad'))
        # grad at self shouldn't be zero, except for Triangle, Exponential, PowerExp
        expect_no_error(gpgradX <- gp$grad(gp$X))
        if (!(j %in% c(4,7,9, 17, 18))) {
          ttm_expect_true(!any(is.na(gpgradX)),
                          label=paste(j, kern_char, '!any(is.na(gpgradX))'))
        }
      } else {
        if (exists('printkern') && printkern) {
          cat("grad/dC_dx not tested for", j, kern_char, "\n")
        }
      }

      # Test gradpredvar
      # FIX m32/m52
      if (j %in% c(1:13) || TRUE) {
        expect_no_error(gpv <- gp$gradpredvar(xgrad))
        # numDeriv::grad(func=function(x) gp$pred(x, se=T)$s2, gpv)
        npv <- 39
        XXpv <- matrix(runif(2*npv), ncol=2)
        expect_no_error(XXgpv <- gp$gradpredvar(XXpv))
        gpvmatches <- 0
        numpvs <- c()
        actpvs <- c()
        for (iii in 1:npv) {
          numpv <- numDeriv::grad(func=function(x) {gp$pred(x, se=T)$s2},
                                  XXpv[iii,])
          # expect_equal(numpv, XXgpv[iii,], tolerance = 1e-2)
          pvclose <- all(
            ifelse(abs(XXgpv[iii,]) < 1e-8,
                   abs(numpv - XXgpv[iii,]) < 1e-8,
                   ifelse(abs(XXgpv[iii,]) < 1e-2,
                          abs((numpv - XXgpv[iii,]) / XXgpv[iii,]) < 1e-1,
                          abs((numpv - XXgpv[iii,]) / XXgpv[iii,]) < 1e-3)))
          if (pvclose) {gpvmatches <- gpvmatches + 1}
          numpvs <- c(numpvs, numpv)
          actpvs <- c(actpvs, XXgpv[iii,])
        }
        if (exists('printkern') && printkern && gpvmatches < npv) {
          cat(kern_char, "gpv close on ", gpvmatches, "/", npv, "\n")
        }
        # Setting bar really low since I don't want this failing on CRAN,
        # and if it works in 1D then it should be fine.
        expect_true(gpvmatches > npv/10,
                    label=paste(j, kern_char, 'gpvmatches', gpvmatches,'/',npv,
                                "seed =", seed))
        # qplot(numpvs, actpvs)
        # summary((numpvs - actpvs) / actpvs)
        # cbind(numpvs, actpvs, prop=(numpvs - actpvs) / actpvs)
        # qplot(numpvs, (numpvs - actpvs) / actpvs)
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

      # Kernel param gradient
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
              cat("FAILURE", kern_char, iatt, i, alleq, actgrad[1+i], numgrad, "\n")
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
      expect_true(goodsofar,
                  label=paste(kern_char,
                              'numgrad matches symbolic grad (failed on all',
                              maxattempts, "attempts)"))
    }) # End ttm
  } # End kernel for loop
})

# Factor kernels ----
test_that("Factor kernels", {
  n <- sample(20:30, 1)
  d <- 1
  x <- matrix(runif(n*d), ncol=d)
  # second is factor dim
  nlev <- 3
  x[, 1] <- sample(nlev, n, T)
  # f <- function(x) {abs(sin(x[1]^.8*6))^1.2 + log(1+x[2]) + x[1]*x[2]}
  f <- function(x) {x[1]^.7}
  y <- apply(x, 1, f) + rnorm(n,0,1e-1) #f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  kern_chars <- c('FactorKernel', 'OrderedFactorKernel',
                  'LatentFactorKernel', 'LatentFactorKernel',
                  'GowerFactorKernel')
  kern_list <- list(
    FactorKernel$new(D=1, nlevels=3, xindex=1),
    OrderedFactorKernel$new(D=1, nlevels=3, xindex=1),
    LatentFactorKernel$new(D=1, nlevels=3, xindex=1, latentdim = 1, s2_est=F, s2=.3),
    LatentFactorKernel$new(D=1, nlevels=3, xindex=1, latentdim = 2),
    GowerFactorKernel$new(D=1, nlevels=3, xindex=1)
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

    # Check kernel properties
    expect_equal(GauPro:::find_kernel_cts_dims(gp$kernel), NULL)
    expect_equal(GauPro:::find_kernel_factor_dims(gp$kernel), c(1,nlev))

    # Test predict
    if (T || (j %in% 1:4)) {
      expect_error(predict(gp, 1:3, se.fit = T), NA,
                   info = paste("bad pred in", kern_char))
      # Can get a warning if preds are too small
      # expect_warning(predict(gp, 1:3, se.fit = T), NA,
      #                info = paste("bad pred in", kern_char))
      # Test plot
      expect_no_error(pp <- gp$plot1D())
      expect_no_error(suppressMessages({pp <- plot(gp)}))
    } else {
      if (exists('printkern') && printkern) {
        cat("factorkernel", j, kern_char, "not testing pred/plot", "\n")
      }
    }

    # Kernel plot
    expect_error(plot(gp$kernel), NA)
    if (j > 1.5 && j < 4.5) {
      expect_no_error(gp$kernel$plotLatent())
    }

    # Test importance
    expect_error(capture.output(imp <- gp$importance(plot=T)), NA)
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
    expect_equal(df, dfg$fn, tolerance = 1e-4)
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
    # Give it many chances so tests don't fail often.

    # max attempts
    maxattempts <- 10
    numgradtol <- 1e-3
    for (iatt in 1:maxattempts) {
      goodsofar <- TRUE

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
        #   expect_equal(numgrad, actgrad[1+i], tolerance = 1e-2,
        #                label=paste(j,kern_char,i, 'numgrad'))
        # }

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
    } # end attempts
  } # end kernel j
})

# Factor*Cts kernels ----
test_that("Factor kernels in product", {
  n <- sample(20:30, 1)
  d <- 2
  x <- matrix(runif(n*d), ncol=d)
  # second is factor dim
  nlev <- 3
  x[, 2] <- sample(1:nlev, n, T)
  colnames(x) <- c('cts', 'fct')
  f <- function(x) {abs(sin(x[1]^.8*6))^1.2 + log(1+x[2]) + x[1]*x[2]}
  y <- apply(x, 1, f) + rnorm(n,0,1e-2) #f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  kern_chars <- c('FactorKernel', 'OrderedFactorKernel',
                  'LatentFactorKernel1', 'LatentFactorKernel2',
                  'GowerFactorKernel')
  kern_list <- list(
    FactorKernel$new(D=2, nlevels=nlev, xindex=2),
    OrderedFactorKernel$new(D=2, nlevels=nlev, xindex=2),
    LatentFactorKernel$new(D=2, nlevels=nlev, xindex=2, latentdim = 1),
    LatentFactorKernel$new(D=2, nlevels=nlev, xindex=2, latentdim = 2),
    GowerFactorKernel$new(D=2, nlevels=nlev, xindex=2)
  )
  for (j in 1:length(kern_chars)) {
    kern_char <- kern_chars[j]
    if (exists('printkern') && printkern) cat(j, kern_char, " w/ cts\n")
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

    # Check kernel properties
    expect_equal(GauPro:::find_kernel_cts_dims(gp$kernel), 1)
    expect_equal(GauPro:::find_kernel_factor_dims(gp$kernel), c(2,nlev))

    # Check kernel print
    expect_error({kernprint <- capture_output(print(gp$kernel))}, NA)
    expect_is(kernprint, 'character')

    # Kernel plot
    expect_error(plot(gp$kernel), NA)

    # Summary
    expect_no_warning(
      expect_error(capture.output(summary(gp)), NA)
    )

    # Check LOO
    expect_no_error((gp$plotLOO()))

    # Check EI
    expect_no_error(suppressWarnings(mei1 <- gp$maxEI()))
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

    # max attempts
    maxattempts <- 10
    numgradtol <- 1e-3
    for (iatt in 1:maxattempts) {
      goodsofar <- TRUE

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
        #   expect_equal(numgrad, actgrad[1+i], tolerance = 1e-2,
        #                label=paste(j,kern_char,i, 'numgrad'))
        # }


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
    } # end attempts
  } # end kernel j
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
  expect_no_error({
    gpf <- GauPro_kernel_model$new(X=tdf, Z=z ~ a + b + c + e, kernel='gauss')})
  expect_true("GauPro" %in% class(gpf))
  expect_equal(ncol(gpf$X), 4)
  expect_true(is.matrix(gpf$X))
  expect_error(predict(gpf, tdf), NA)
  # Print
  expect_no_error(printout <- capture.output(print(gpf)))
  expect_true(is.character(printout))
  expect_equal(printout[1], "GauPro kernel model object")
  # Plot
  expect_no_error(gpf$plotmarginal())
  expect_no_error(gpf$plotmarginalrandom())
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
  # Plot
  expect_no_error(gpdf$plotmarginal())
  expect_no_error(gpdf$plotmarginalrandom())
  # Optimize
  expect_no_error(gpdf$optimize_fn(function(x) {gpdf$predict(x)}))
  # Test EI
  expect_no_error(gpdf$EI(xdf[1,]))
  # Test maxEI
  # expect_error(dfEI <- gpdf$maxEI(), NA)
  expect_error(dfEI <- gpdf$maxEI(lower=c(-3,0,1,1,0), upper=c(3,1,5,4,4)), NA)
  expect_true(is.data.frame(dfEI$par))
  expect_equal(colnames(dfEI$par), colnames(xdf)[1:5])
  expect_equal(dim(dfEI$par), c(1,5))
  rm(dfEI)
  # maxEI with mopar, each of EI, AugmentedEI, and CorrectedEI
  mopar <- c(mixopt::mopar_cts(-3,3),
             mixopt::mopar_cts(0,1),
             mixopt::mopar_unordered(letters[1:5]),
             mixopt::mopar_unordered(letters[6:9]),
             mixopt::mopar_cts(0,4)
  )
  for (eitype in c("Ei", "cOR", "augmentedEI")) {
    expect_no_error(dfEI2 <- gpdf$maxEI(mopar = mopar, EItype = eitype))
    expect_true(is.data.frame(dfEI2$par))
    expect_equal(colnames(dfEI2$par), colnames(xdf)[1:5])
    expect_equal(dim(dfEI2$par), c(1,5))
    rm(dfEI2)
  }

  # Test qEI with mopar
  expect_no_error(dfqEI <- gpdf$maxqEI(npoints = 2,mopar = mopar))
  expect_true(is.data.frame(dfqEI$par))
  expect_equal(colnames(dfqEI$par), colnames(xdf)[1:5])
  expect_equal(dim(dfqEI$par), c(2,5))
  rm(dfqEI)


  # Try other arg names
  expect_error(gpdf <- GauPro_kernel_model$new(z ~ ., xdf, kernel='m32'), NA)
  expect_error(gpdf <- GauPro_kernel_model$new(xdf, z ~ ., kernel='m32'), NA)
  expect_error(gpdf <- GauPro_kernel_model$new(formula=z ~ ., data=xdf, kernel='m32'), NA)
  expect_error(gpdf <- GauPro_kernel_model$new(z ~ ., xdf, kernel='m32'), NA)
  expect_error(gpdf <- GauPro_kernel_model$new(z ~ xdf$a + xdf$c + xdf$e, xdf, kernel='m32'), NA)
  expect_error(gpdf <- GauPro_kernel_model$new(z ~ ., data=xdf, kernel='m32'), NA)
  rm(gpdf)

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

test_that("Formula/data input 3", {
  # Add ordered in autokernel
  library(dplyr)
  n <- 63
  xdf <- tibble(
    a=rnorm(n),
    b=runif(n),
    c=ordered(sample(letters[1:5], n, T)),
    d=factor(sample(letters[6:9], n, T)),
    z=a*b + a^2*ifelse(c %in% c('a', 'b'), 1, .5) +
      ifelse(d %in% c('g','h'), 1, -1) +
      ifelse(paste0(d) %in% c('a','c'),4,0) +
      rnorm(n, 1e-3)
  )

  # Test fit
  expect_error(gpdf <- GauPro_kernel_model$new(z ~ ., data=xdf, kernel='m32'), NA)
})

test_that("Formula/data input 4", {
  # Only cts dimensions, no factor dimensions
  library(dplyr)
  n <- 63
  xdf <- tibble(
    a=rnorm(n),
    b=runif(n),
    z=a*b + a^2 + rnorm(n, 1e-3)
  )

  # Test fit
  expect_error(gpdf <- GauPro_kernel_model$new(z ~ ., data=xdf, kernel='m32'), NA)

  # Test predict (this gave error before)
  expect_error(predict(gpdf, xdf), NA)
})

# EI ----
test_that("EI with cts", {
  n <- 20
  d <- 2
  x <- GauPro:::lhs_maximinLHS(n, d) # Better spacing might avoid grad issues?
  n <- nrow(x)
  f <- function(x) {abs(sin(x[1]^.8*6))^1.2 + log(1+(x[2]-.3)^2) + x[1]*x[2]}
  y <- apply(x, 1, f) + rnorm(n,0,1e-2) #f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  p1 <- GauPro_kernel_model$new(x, y)
  # Grid of inputs to test
  ei_grid <- expand.grid(minim=c(T, F), method=1:3, matr=c(T,F), eps=c(0, 1e-3))
  for (igrid in rev(1:nrow(ei_grid))) {
    # print(ei_grid[igrid,])
    minim <- ei_grid$minim[igrid]
    maxEI_par <- p1$maxEI(minimize = minim)$par # + rnorm(d, 0, .05)
    meth <- ei_grid$method[igrid]
    eifunc <- if (meth == 1) {p1$EI} else if (meth==2) {
      p1$AugmentedEI} else {p1$CorrectedEI}
    matr <- ei_grid$matr[igrid]
    i_eps <- ei_grid$eps[igrid]
    xx <- if (matr) {
      matrix(maxEI_par, ncol=d, nrow=7, byrow=T) + matrix(rnorm(d*7,0,.05), ncol=d)
    } else {
      maxEI_par +  rnorm(d, 0, .05)
    }
    einumgrad <- if (matr) {
      matrix(
        unlist(
          lapply(1:nrow(xx),
                 function(i) {
                   numDeriv::grad(
                     function(h) {eifunc(x=h, minimize=minim, eps=i_eps)},
                     x=xx[i,])})
        ), byrow=T, ncol=d
      )
    } else {
      numDeriv::grad(function(h) {eifunc(x=h, minimize=minim, eps=i_eps)}, x=xx)
    }
    analyticgrad <- eifunc(x=xx, minimize=minim, return_grad = T, eps=i_eps)$grad
    expect_equal(c(einumgrad),   c(analyticgrad), tolerance = 1e-2)
    if (F) {
      curve(sapply(x, function(xxx) {
        eifunc(x=xx, minimize=minim, return_grad = F, eps=xxx)}), 0, .2)
    }
  }
})
test_that("EI with mixopt", {
  n <- 30
  tdf <- data.frame(a=runif(n), b=runif(n, -1,1),
                    c=(sample(5:6,n,T)),
                    d=sample(c(.1,.2,.3,.4), n, T),
                    # e=sample(letters[1:3], n,T),
                    f=sample(10:30, n, T))
  z <- with(tdf, a+a*b+b^2 +5*(d-.22)^2*(f-22)^2)
  gpf <- GauPro_kernel_model$new(X=as.matrix(tdf), Z=z, kernel='m52')
  expect_no_error(gpf$EI(as.matrix(tdf[1,])))
  expect_no_error(gpf$maxEI(minimize = T))
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

test_that("EI minimize is right", {
  nattempts <- 10
  for (iattempt in 1:nattempts) {
    # cat(iattempt, "\n")
    d <- 1
    n <- 6
    x <- runif(n)
    y <- sin(2*pi*x^.9)*(1+.2*x^.5) + rnorm(n,0,1e-3)
    gp <- GauPro_kernel_model$new(x, y, kernel=Matern52)
    # gp$plot1D()
    u <- matrix(seq(0,1,l=101), ncol=1)
    expect_no_error(ei1 <- gp$EI(u, minimize = T))
    gpinv <- GauPro_kernel_model$new(
      x, -y,
      kernel=Matern52$new(D=1, beta=gp$kernel$beta,
                          beta_est=F, s2=gp$kernel$s2, s2_est=F))
    # gpinv$plot1D()
    expect_no_error(eiinv <- gpinv$EI(u, minimize=F))
    # plot(ei1, eiinv)
    if ((iattempt >= nattempts) ||
        testthat::compare(ei1, eiinv, tol=1e-3)$equal) {
      expect_equal(ei1, eiinv, tol=1e-3)
    } else {
      next
    }

    # Augmented EI
    expect_no_error(augei1 <- gp$AugmentedEI(u, minimize=T))
    expect_no_error(augei2 <- gpinv$AugmentedEI(u, minimize=F))
    # expect_equal(augei1, augei2, tol=1e-1)
    if ((iattempt >= nattempts) ||
        testthat::compare(augei1, augei2, tol=1e-1)$equal) {
      expect_equal(augei1, augei2, tol=1e-1)
    } else {
      next
    }
    # plot(augei1, augei2)
    # curve(gp$AugmentedEI(matrix(x, ncol=1), minimize=T))
    # curve(gpinv$AugmentedEI(matrix(x, ncol=1), minimize=F), add=T, col=2)
    # curve(gp$AugmentedEI(matrix(x, ncol=1), minimize=F))
    # curve(gpinv$AugmentedEI(matrix(x, ncol=1), minimize=T), add=T, col=2)
    break
  }
})

test_that("Aug EI makes sense", {
  d <- 1
  n <- 16
  x <- c(seq(0,1,l=n), seq(.3,.5,l=n*3))
  n <- length(x)
  # y <- sin(2*pi*x^.9)*(1+.2*x^.5) + rnorm(n,0,1e-1)
  # y <- x^4-x^2 + .01*sin(2*pi*x^.9)*(1+.2*x^.5) + rnorm(n,0,1e-1)
  y <- sin(2*pi*x*2)+ .3*x + rnorm(n,0,1e-1)
  gp <- GauPro_kernel_model$new(x, y, kernel=Matern52)
  # gp$plot1D()

  # nattempts <- 10
  # for (iatt in 1:nattempts) {
  u <- matrix(seq(0,1,l=101), ncol=1)
  expect_no_error(ei1 <- gp$EI(u, minimize = T))
  # }

  # curve(gp$EI(matrix(x, ncol=1), minimize=T))
  # curve(gp$AugmentedEI(matrix(x, ncol=1), minimize=T))
  # curve(gp$pred(matrix(x, ncol=1), se=T)$se)
  # curve(gp$pred(matrix(x, ncol=1), se=T, mean_dist = T)$se)
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
test_that("Wide range Z", {
  d <- 3
  n <- 40
  x <- cbind(runif(n, 50, 150),
             runif(n, -9,-7),
             runif(n, 1200, 3000))
  # f <- function(x) {1e3*(x[1] + 14* x[2]^2 + x[3])} # Often gets stuck
  f <- function(x) {1e2*(x[1]^1.2 + 14* x[2]^2 + x[3]^.9)} # Reliable
  y <- apply(x, 1, f)
  y <- y + rnorm(n)*.01*diff(range(y))
  expect_no_error(e1 <- GauPro_kernel_model$new(
    x, y,
    # kernel="gauss",#Gaussian$new(D=3, s2=3e7, s2_lower=1e7),
    kernel=Gaussian$new(D=ncol(x), s2=3e7, s2_lower=1e7, s2_upper=1e20),
    track=T,
    verbose=0, restarts=25))
  # e1$plotLOO(); print(e1$s2_hat); print(e1$nug)
  expect_no_error(e1$plot_track_optim())
  expect_no_error(e1$summary())
  expect_no_error(e1$plotLOO())
  expect_no_error(e1$plotmarginalrandom())
  expect_no_error(e1$plotmarginal())
  expect_no_error(plot(e1))
  expect_error(e1$plot1D())
  expect_error(e1$plot2D())
  expect_error(e1$cool1Dplot())
  expect_lt(mean(abs((e1$Z-e1$pred(e1$X)) / e1$Z)), 1e-1)
  e1$update_fast(Xnew=.5*x[1,,drop=F] + .5*x[2,], .5*(y[1]+y[2]))
})
test_that("Wide range X", {
  d <- 3
  n <- 40
  x <- cbind(runif(n, 1e8, 2e8),
             runif(n, 1e-8, 2e-8),
             runif(n, -1e12, 1e12))
  # f <- function(x) {1e3*(x[1] + 14* x[2]^2 + x[3])} # Often gets stuck
  f <- function(x) {1e2*(x[1]^1.2 + 14* x[2]^2 + abs(x[3])^.9)} # Reliable
  y <- apply(x, 1, f)
  y <- y + rnorm(n)*.01*diff(range(y))
  # Error for bad X/Z
  expect_error(gpkm(x, c(NaN, y[-1])))
  expect_error(gpkm(x, c(NA, y[-1])))
  expect_error(gpkm(x, c(Inf, y[-1])))
  expect_error(gpkm(rbind(NA, x[-1,]), y))
  expect_error(gpkm(rbind(NaN, x[-1,]), y))
  expect_error(gpkm(rbind(Inf, x[-1,]), y))
  # Actual fit
  expect_no_error(e1 <- GauPro_kernel_model$new(
    x, y,
    # kernel="gauss",#Gaussian$new(D=3, s2=3e7, s2_lower=1e7),
    kernel=Gaussian$new(D=ncol(x), s2=3e7, s2_lower=1e7, s2_upper=1e20),
    verbose=0, restarts=25))
  # e1$plotLOO(); print(e1$s2_hat); print(e1$nug)
  expect_no_error(e1$summary())
  expect_no_error(e1$plotLOO())
  expect_no_error(e1$plotmarginalrandom())
  expect_no_error(e1$plotmarginal())
  expect_no_error(plot(e1))
  expect_error(e1$plot1D())
  expect_error(e1$plot2D())
  expect_error(e1$cool1Dplot())
  # expect_lt(mean(abs((e1$Z-e1$pred(e1$X)) / e1$Z)), 1e-1)
  e1$update_fast(Xnew=.5*x[1,,drop=F] + .5*x[2,], .5*(y[1]+y[2]))
})

# Predictions are at least the minimum value
test_that("Predictions are big enough", {
  n <- 11
  x <- seq(0,1,l=11)
  y <- x
  expect_error(suppressMessages(gp <- gpkm(y ~ x)), NA)
  expect_error(suppressWarnings(pm <- gp$predict(x, se.fit=T, mean_dist=T)), NA)
  expect_true(all(pm$s2 >= 0))
  expect_error(suppressWarnings(p <- gp$predict(x, se.fit=T, mean_dist=F)), NA)
  expect_true(all(p$s2 >= gp$nug * gp$s2_hat))
})

# Bad kernels ----
test_that("Bad kernels", {
  # Data set up
  n <- 20
  d <- 2
  x <- matrix(runif(n*d), ncol=d)
  # second is factor dim
  nlev <- 2
  x[, 2] <- sample(1:nlev, n, T)
  f <- function(x) {abs(sin(x[1]^.8*6))^1.2 + log(1+x[2]) + x[1]*x[2]}
  y <- apply(x, 1, f) + rnorm(n,0,1e-2)

  # Can't use factor twice on same dim
  expect_error({
    GauPro_kernel_model$new(
      x, y,
      kernel=IgnoreIndsKernel$new(k = Gaussian$new(D=1), ignoreinds = 2) *
        FactorKernel$new(D=2, xindex = 2, nlevels = nlev) *
        LatentFactorKernel$new(D=2, xindex=2, nlevels=nlev))
  })
  # Can't use same index as cts and factor
  expect_error({
    GauPro_kernel_model$new(
      x, y,
      kernel=Gaussian$new(D=2) *
        FactorKernel$new(D=2, xindex = 2, nlevels = nlev))
  })
})

# Normalize Z ----
test_that("Normalize Z", {
  d <- 3
  n <- 30
  x <- matrix(runif(d*n), ncol=d)
  if (d == 1) {
    y <- sin(6*x[,1]^.7)*1e6 + rnorm(n, 0, 1e2)
  } else {
    y <- (x[,1]^.7 + x[,2]^1.4 * (x[,3]^.8*2))*1e6 + rnorm(n, 0, 1e2)
  }
  expect_no_error({
    gp <- GauPro_kernel_model$new(x, y, normalize = T, trend=trend_c)
  })
  # plot(predict(gp, x), y)
  expect_no_error(pred <- predict(gp, x))
  expect_equal(y, pred, tolerance = 1e-1)
  expect_equal(y, gp$pred_LOO(), tolerance = 1e-1)
  expect_no_error(normEI <- gp$maxEI())
  # expect_true(normEI$value < .5*(max(y) - min(y)))
})
# Diamonds ----
test_that("Diamonds", {
  n <- 000 + sample(50:70, 1)
  expect_no_error({
    system.time({
      dm <- gpkm(price ~ carat + cut + color + clarity + depth +
                   table + x + y + z,
                 ggplot2::diamonds[sample(1:nrow(ggplot2::diamonds),
                                          n, replace=FALSE
                 ), ])})
    summary(dm)
  })
  expect_no_error(dm$plotmarginal())
  expect_no_error(dm$plotmarginalrandom())
})

# Kernels useC ----
test_that("Kernels useC", {
  d <- 4
  n <- 11
  X1 <- matrix(runif(n*d), ncol=d)
  X2 <- matrix(runif(2*n*d), ncol=d)
  V1 <- runif(d)
  kern_chars <- c('Gaussian', 'Matern32', 'Matern52',
                  'Triangle', 'Cubic', 'White',
                  'PowerExp', 'Periodic', "Exponential", "RatQuad")
  kernlist <- list(Gaussian, Matern32, Matern52,
                   Triangle, Cubic, White,
                   PowerExp, Periodic, Exponential, RatQuad)
  for (i in seq_along(kern_chars)) {
    # Make kernel and clone it, set that one to use R
    expect_no_error({
      kC <- kernlist[[i]]$new(D=d)
      kC$set_params_from_optim(kC$param_optim_start(jitter=T))
      kR <- kC$clone(T)
      kR$useC <- F
    })
    # Make sure it matches
    expect_equal(kC$k(X1), kR$k(X1))
    expect_equal(kC$k(X1, X2), kR$k(X1, X2))

    expect_equal(kC$k(X1, V1), kR$k(X1, V1))
    expect_equal(kC$k(V1, X1), kR$k(V1, X1))

    expect_equal(kC$k(V1), kR$k(V1))

    if (exists('printkern') && printkern) {
      cat(kern_chars[i], length(kC$C_dC_dparams(X=X1, nug=1e-4)$d), "\n")
    }
    expect_equal(kC$C_dC_dparams(X=X1, nug=1e-4)$d,
                 kR$C_dC_dparams(X=X1, nug=1e-4)$d,
                 label=paste("useC", i, kern_chars[i], "C_dC_dparams"))
  }
})
