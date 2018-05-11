test_that("kernel_model works", {
  # No longer using this old Gaussian kernel since it wasn't on log scale.
  # set.seed(0)
  # n <- 20
  # x <- matrix(seq(0,1,length.out = n), ncol=1)
  # f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
  # y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  # gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian$new(1), parallel=FALSE, verbose=10, nug.est=T)
  # # gp$cool1Dplot()
  # # numDeriv::grad(func = gp$deviance, x=c(5,1))
  # # gp$deviance_grad(params = c(5,1), nug.update=F)
  #
  # expect_equal(gp$kernel$theta, 4.826752, tolerance=.1)
  # expect_equal(gp$kernel$s2, 1.963462, tolerance=.1)
  # expect_equal(gp$nug, 0.000435802, tolerance=.01)
  # expect_equal(
  #   # numDeriv::grad(func = function(x) gp$deviance(params=x[1:2], nuglog=x[3]), x=c(5,1, -4)),
  #   c(-6.641381,  -98.038047, -213.498400),
  #   gp$deviance_grad(params = c(5,1), nug.update=T, nuglog=-4),
  #   tol=.01
  #   )
})
test_that("kernel_Gaussian works", {

  set.seed(0)
  n <- 20
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
  y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian$new(1), parallel=FALSE, verbose=10, nug.est=T)
  # gp$cool1Dplot()
  # numDeriv::grad(func = gp$deviance, x=c(5,1))
  # gp$deviance_grad(params = c(5,1), nug.update=F)

  expect_equal(gp$kernel$beta, 1.185298, tolerance=.1)
  expect_equal(gp$kernel$s2, 0.3809001, tolerance=.1)
  expect_equal(gp$nug, 0.001037888, tolerance=.0001)
  expect_equal(
    # numDeriv::grad(func = function(x) gp$deviance(params=x[1:2], nuglog=x[3]), x=c(.2,1.1, -4.5)),
    c(-272.12230, -102.95356,  -63.01746),
    gp$deviance_grad(params = c(.2,1.1), nug.update=T, nuglog=-4.5, trend_update=FALSE),
    tol=.001
  )
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:2], nuglog=x[3])}, x=c(2.5,-.2, -5)),
    c(30.89216, 26.80619,  0.0006299453),
    gp$deviance_grad(params = c(2.5,-.2), nug.update=T, nuglog = -5, trend_update=FALSE),
    tol=.001
  )
})
# test_that("kernel_Gaussian_l works", {
#   set.seed(0)
#   n <- 20
#   x <- matrix(seq(0,1,length.out = n), ncol=1)
#   f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
#   y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
#   gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian_l$new(1), parallel=FALSE, verbose=10, nug.est=T)
#
#   expect_equal(gp$kernel$l, 0.1806493, tolerance=.01)
#   expect_equal(gp$kernel$s2, 0.3809001, tolerance=.01)
#   expect_equal(gp$nug, 0.001037888, tolerance=.0001)
#   expect_equal(
#     # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:2], nuglog=x[3])}, x=c(-.7,.227, -3.66)),
#     c(-20.509641,  11.177885,  -1.592343),
#     gp$deviance_grad(params = c(-.7,.227), nug.update=T, nuglog=-3.66),
#     tol=.001
#   )
#   expect_equal(
#     # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:2], nuglog=x[3])}, x=c(2.5,-.2, -5)),
#     c(3847404.2, -2317793.3, -1894093.7),
#     gp$deviance_grad(params = c(2.5,-.2), nug.update=T, nuglog = -5),
#     tol=100
#   )
# })
test_that("kernel_Exponential works", {
  set.seed(0)
  n <- 20
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
  y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Exponential$new(1), parallel=FALSE, verbose=10, nug.est=T)

  expect_equal(gp$kernel$beta, 0.430569, tolerance=.01)
  expect_equal(gp$kernel$s2, 0.3158615, tolerance=.01)
  expect_equal(gp$nug, 2.846939e-11, tolerance=.0001)
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:2], nuglog=x[3])}, x=c(-.7,.227, -3.66)),
    c(6.4639755,  16.3037264,  0.3593728),
    gp$deviance_grad(params = c(-.7,.227), nug.update=T, nuglog=-3.66, trend_update=FALSE),
    tol=.001
  )
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:2], nuglog=x[3])}, x=c(2.5,-.2, -5)),
    c(13.011660214, 28.946608363, 0.000535044),
    gp$deviance_grad(params = c(2.5,-.2), nug.update=T, nuglog = -5, trend_update=FALSE),
    tol=.001
  )
})
test_that("kernel_Matern32 works", {
  set.seed(0)
  n <- 20
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
  y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Matern32$new(1), parallel=FALSE, verbose=10, nug.est=T, nug.min=0)

  expect_equal(gp$kernel$beta, 0.4609827, tolerance=.01)
  expect_equal(gp$kernel$s2, 1.083443, tolerance=.01)
  expect_equal(log(gp$nug,10), -11.37841, tolerance=.1)
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:2], nuglog=x[3])}, x=c(-.7,.227, -3.66)),
    c(-877.52290, -643.49553, -30.76981),
    gp$deviance_grad(params = c(-.7,.227), nug.update=T, nuglog=-3.66, trend_update=FALSE),
    tol=.01
  )
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:2], nuglog=x[3])}, x=c(2.5,-.2, -5)),
    c(2.393571e+01, 3.066583e+01, 7.917006e-04),
    gp$deviance_grad(params = c(2.5,-.2), nug.update=T, nuglog = -5, trend_update=FALSE),
    tol=.001
  )
})
test_that("kernel_Matern52 works", {
  set.seed(0)
  n <- 20
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
  y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Matern52$new(1), parallel=FALSE, verbose=10, nug.est=T)

  expect_equal(gp$kernel$beta, 0.7997398, tolerance=.01)
  expect_equal(gp$kernel$s2, 0.9536505, tolerance=.01)
  expect_equal(log(gp$nug,10), -3.490174, tolerance=.1)
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:2], nuglog=x[3])}, x=c(-.7,.227, -3.66)),
    c(-9005.369, -5128.656, -1037.114),
    gp$deviance_grad(params = c(-.7,.227), nug.update=T, nuglog=-3.66, trend_update=FALSE),
    tol=.1
  )
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:2], nuglog=x[3])}, x=c(2.5,-.2, -5)),
    c(30.80262739, 31.13268023, 0.00100013),
    gp$deviance_grad(params = c(2.5,-.2), nug.update=T, nuglog = -5, trend_update=FALSE),
    tol=.001
  )
})
test_that("kernel_sum works", {
  set.seed(0)
  n <- 20
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
  y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian$new(.7)+Matern32$new(1.2), parallel=FALSE, verbose=10, nug.est=T)

  expect_equal(gp$kernel$k1$beta, 0.7717338, tolerance=.01)
  expect_equal(gp$kernel$k1$s2, 0.93595, tolerance=.01)
  expect_equal(gp$kernel$k2$beta, 1.70821, tolerance=.01)
  expect_equal(gp$kernel$k2$s2, 0.004285955, tolerance=.01)
  expect_equal(log(gp$nug,10), -3.531732, tolerance=.1)
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:4], nuglog=x[5])}, x=c(-.7,.227,1.1,.3, -3.66)),
    c(0.219613,  2.603450, 45.827544, 38.221117,  1.613496),
    gp$deviance_grad(params = c(-.7,.227,1.1,.3), nug.update=T, nuglog=-3.66, trend_update=FALSE),
    tol=.001
  )
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:4], nuglog=x[5])}, x=c(2.5,-.2, -.9,-.2, -5)),
    c(30.065840109, 28.008131903, -1.052814855,  0.993413933,  0.001268154),
    gp$deviance_grad(params = c(2.5,-.2,-.9,-.2), nug.update=T, nuglog = -5, trend_update=FALSE),
    tol=.001
  )

  # Again with different kernels
  set.seed(0)
  n <- 20
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
  y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Exponential$new(.7)+Matern52$new(1.2), parallel=FALSE, verbose=10, nug.est=T)

  expect_equal(gp$kernel$k1$beta, 5.186118, tolerance=.5)
  expect_equal(gp$kernel$k1$s2, 0.0002988821, tolerance=.01)
  expect_equal(gp$kernel$k2$beta, 0.7997681, tolerance=.01)
  expect_equal(gp$kernel$k2$s2, 0.9535706, tolerance=.01)
  expect_equal(log(gp$nug,10), -4.997353, tolerance=.1)
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:4], nuglog=x[5])}, x=c(-.7,.227,1.1,.3, -3.66)),
    c(13.4913194, 29.1705405, 15.2660073, 12.2975015, 0.7419428),
    gp$deviance_grad(params = c(-.7,.227,1.1,.3), nug.update=T, nuglog=-3.66, trend_update=FALSE),
    tol=.001
  )
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:4], nuglog=x[5])}, x=c(2.5,-.2, -.9,-.2, -5)),
    c(13.029202435, 28.573028878, -0.604253831,  1.499179095,  0.001062688),
    gp$deviance_grad(params = c(2.5,-.2,-.9,-.2), nug.update=T, nuglog = -5, trend_update=FALSE),
    tol=.001
  )
})
test_that("kernel_product works", {
  set.seed(0)
  n <- 20
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
  y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian$new(.7,.01)*Matern32$new(1.2,1.8), parallel=FALSE, verbose=10, nug.est=T)

  expect_equal(gp$kernel$k1$beta, 0.7671658, tolerance=.01)
  expect_equal(gp$kernel$k1$s2, 0.07390575, tolerance=.01)
  expect_equal(gp$kernel$k2$beta, -0.0775485, tolerance=.01)
  expect_equal(gp$kernel$k2$s2, 10.60559, tolerance=.01)
  expect_equal(log(gp$nug,10), -3.344504, tolerance=.1)
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:4], nuglog=x[5])}, x=c(-.7,.227,1.1,.3, -3.66)),
    c(0.2117085, 42.9840666, 47.8857571, 42.9840666,  0.9225288),
    gp$deviance_grad(params = c(-.7,.227,1.1,.3), nug.update=T, nuglog=-3.66, trend_update=FALSE),
    tol=.001
  )
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:4], nuglog=x[5])}, x=c(2.5,-.2, -.9,-.2, -5)),
    c(36.06449032, 15.57771736,  0.02086093, 15.57771736,  0.00056776),
    gp$deviance_grad(params = c(2.5,-.2,-.9,-.2), nug.update=T, nuglog = -5, trend_update=FALSE),
    tol=.001
  )

  # Again with different kernels
  set.seed(0)
  n <- 20
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
  y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Exponential$new(.7)*Matern52$new(1.2), parallel=FALSE, verbose=10, nug.est=T)

  expect_equal(gp$kernel$k1$beta, -20.25548, tolerance=5) # Tiny theta, give it large tolerance
  expect_equal(gp$kernel$k1$s2, 0.976513, tolerance=.01)
  expect_equal(gp$kernel$k2$beta, 0.7997657, tolerance=.01)
  expect_equal(gp$kernel$k2$s2, 0.976513, tolerance=.01)
  expect_equal(log(gp$nug,10), -3.49013, tolerance=.1)
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:4], nuglog=x[5])}, x=c(-.7,.227,1.1,.3, -3.66)),
    c(14.2077686, 43.3416492, 14.9012499, 43.3416491,  0.3427752),
    gp$deviance_grad(params = c(-.7,.227,1.1,.3), nug.update=T, nuglog=-3.66, trend_update=FALSE),
    tol=.001
  )
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:4], nuglog=x[5])}, x=c(2.5,-.2, -.9,-.2, -5)),
    c(1.626385e+01, 1.892872e+01, 1.355478e-02, 1.892872e+01, 4.854032e-04),
    gp$deviance_grad(params = c(2.5,-.2,-.9,-.2), nug.update=T, nuglog = -5, trend_update=FALSE),
    tol=.001
  )
})
test_that("kernel_RatQuad works", {
  set.seed(0)
  n <- 20
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
  y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=RatQuad$new(1, alpha=1), parallel=FALSE, verbose=10, nug.est=T)

  expect_equal(gp$kernel$beta, 1.10732, tolerance=.01)
  expect_equal(gp$kernel$alpha, 12.37519, tolerance=.01)
  expect_equal(gp$kernel$logs2, -0.365198, tolerance=.01)
  expect_equal(log(gp$nug,10), -3.0309, tolerance=.1)
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:3], nuglog=x[4])}, x=c(-.7,.227, .65, -3.66)),
    c(-6659.246, 2182.775, -3746.691, -1451.918),
    gp$deviance_grad(params = c(-.7,.227, .65), nug.update=T, nuglog=-3.66, trend_update=FALSE),
    tol=.1
  )
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:3], nuglog=x[4])}, x=c(2.5, 1.2, -.2, -5)),
    c(3.000658e+01, 7.499106e-01, 2.720578e+01, 6.428928e-04),
    gp$deviance_grad(params = c(2.5, 1.2, -.2), nug.update=T, nuglog = -5, trend_update=FALSE),
    tol=.01
  )
})
test_that("kernel_Periodic works", {
  set.seed(0)
  n <- 20
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
  y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Periodic$new(p=1, alpha=1), parallel=FALSE, verbose=10, nug.est=T)

  expect_equal(gp$kernel$logp, 0.04179651, tolerance=.01)
  expect_equal(gp$kernel$logalpha, 1.045002, tolerance=.01)
  expect_equal(gp$kernel$logs2, -0.3807649, tolerance=.01)
  expect_equal(log(gp$nug,10), -3.017538, tolerance=.1)
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:3], nuglog=x[4])}, x=c(-.7,.227, .65, -3.66)),
    c(-2853.328, -1118.391, -10108.507, -9670.486),
    gp$deviance_grad(params = c(-.7,.227, .65), nug.update=T, nuglog=-3.66, trend_update=FALSE),
    tol=.1
  )
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:3], nuglog=x[4])}, x=c(.5, 1.2, -.2, -5)),
    c(-14015.882197, -475.132126, -82.258573, -5.474445),
    gp$deviance_grad(params = c(.5, 1.2, -.2), nug.update=T, nuglog = -5, trend_update=FALSE),
    tol=.01
  )
})
test_that("kernel_PowerExp works", {
  set.seed(0)
  n <- 20
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
  y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=PowerExp$new(1, alpha=1), parallel=FALSE, verbose=10, nug.est=T, restarts=1)

  expect_equal(gp$kernel$beta, 0.783801906330017, tolerance=.01)
  expect_equal(gp$kernel$alpha, 1.98088200644922, tolerance=.01)
  expect_equal(gp$kernel$logs2, -0.0514541683404154, tolerance=.01)
  expect_equal(log(gp$nug,10), -3.565321, tolerance=.1)
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:3], nuglog=x[4])}, x=c(-.7,.227, .65, -3.66)),
    c(27.86004234, -43.18817540,  31.76892239,   0.07266018),
    gp$deviance_grad(params = c(-.7,.227, .65), nug.update=T, nuglog=-3.66, trend_update=FALSE),
    tol=.1
  )
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1:3], nuglog=x[4])}, x=c(2.5, 1.2, -.2, -5)),
    c(0.0600704370, -0.0768150675, 10.8086791014,  0.0001081507),
    gp$deviance_grad(params = c(2.5, 1.2, -.2), nug.update=T, nuglog = -5, trend_update=FALSE),
    tol=.01
  )
})

test_that("kernel_White works", {
  set.seed(0)
  n <- 20
  x <- matrix(seq(0,1,length.out = n), ncol=1)
  f <- Vectorize(function(x) {sin(2*pi*x) + .001*sin(8*pi*x) +rnorm(1,0,.03)})
  y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
  gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=White$new(s2=1), parallel=FALSE, verbose=10, nug.est=F)

  expect_equal(gp$kernel$logs2, -0.3164508, tolerance=.01)
  # expect_equal(log(gp$nug,10), -3.017538, tolerance=.1)
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1])}, x=c(-.7)),
    c(-65.32514),
    gp$deviance_grad(params = c(-.7), nug.update=F, trend_update=FALSE),
    tol=.1
  )
  expect_equal(
    # numDeriv::grad(func = function(x) {gp$deviance(params=x[1])}, x=c(.5)),
    c(39.0243),
    gp$deviance_grad(params = c(.5), nug.update=F, trend_update=FALSE),
    tol=.01
  )
})


