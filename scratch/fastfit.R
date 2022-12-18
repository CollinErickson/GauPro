# Plan was to try to initial fit on subset of data, then go to all data
# for final steps.
# Didn't work, the final step ended up using just as many steps.

fastfit <- function(self) {
  # Save all data
  fullX <- self$X
  fullZ <- self$Z
  # Fit on subset of data
  n1 <- 75
  stopifnot(n1 < nrow(fullX))
  inds1 <- sample(1:nrow(fullX), n1, replace=FALSE)
  X1 <- fullX[inds1, ]
  Z1 <- fullZ[inds1]
  self$X <- X1
  self$Z <- Z1
  self$N <- nrow(self$X)
  # self$update_K_and_estimates()
  self$update_K_and_estimates() # Need to get mu_hat before starting
  system.time(self$update())
  # Now go back to all data
  self$X <- fullX
  self$Z <- fullZ
  self$N <- nrow(self$X)
  # debugonce(self$optim)
  system.time(self$update())
}
