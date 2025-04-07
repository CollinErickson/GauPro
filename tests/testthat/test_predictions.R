test_that("prediction with no categorical predictors", {
  set.seed(1)
  dat <- data.frame(x = runif(20), y = runif(20), z = runif(20))
  gp_kern <- k_Exponential(D = 2)
  gp_fit <- gpkm(z ~ x + y, data = dat, kernel = gp_kern)

  preds <- gp_fit$pred(dat[, 1:2], se.fit = TRUE)

  exp_ptype <-
    structure(
      list(mean = numeric(0), s2 = numeric(0), se = numeric(0)),
      row.names = integer(0),
      class = "data.frame")

  expect_equal(preds[0,], exp_ptype)
  expect_equal(nrow(preds), nrow(dat))

})

test_that("prediction with only categorical predictors", {
  set.seed(1)
  dat <- expand.grid(x = letters[1:2], y = LETTERS[1:2])
  dat$z <- runif(nrow(dat))

  gp_kern <- k_FactorKernel(D = 2, nlevels= 2, xindex = 1:2)
  expect_message(
    gp_fit <- gpkm(z ~ x + y, data = dat, kernel = gp_kern),
    regexp = "All restarts had error, keeping initial"
  )

  preds <- gp_fit$pred(dat[, 1:2], se.fit = TRUE)

  exp_ptype <-
    structure(
      list(mean = numeric(0), s2 = numeric(0), se = numeric(0)),
      row.names = integer(0),
      class = "data.frame")

  expect_equal(preds[0,], exp_ptype)
  expect_equal(nrow(preds), nrow(dat))

})

test_that("prediction with mixed predictor types", {
  # skip("Appears to cause an inf loop of warnings")
  # This runs fine interactively but, during testing, is issues infinite
  # warnings.
  set.seed(1)
  dat <- data.frame(x = runif(20), z = runif(20))
  dat$y <- rep(letters[1:2], 10)


  gp_kern <-  k_Exponential(D = 1) * k_FactorKernel(D = 2, nlevels= 2, xindex = 2)
  gp_fit <- gpkm(z ~ x + y, data = dat, kernel = gp_kern)

  preds <- gp_fit$pred(dat[, c(1, 3)], se.fit = TRUE)

  exp_ptype <-
    structure(
      list(mean = numeric(0), s2 = numeric(0), se = numeric(0)),
      row.names = integer(0),
      class = "data.frame")

  expect_equal(preds[0,], exp_ptype)
  expect_equal(nrow(preds), nrow(dat))

})

