# Augmented EI
# Goal is minimize
# Huang, http://dx.doi.org/10.1007/s10898-005-2454-3
self <- gp1
pred_X <- self$predict(self$X, se.fit = T)
u_X <- -pred_X$mean - pred_X$se
star_star_index <- which.max(u_X)
f <- pred_X$mean[star_star_index]
EI <- function(y, s) {
  z <- (f - y) / s
  (f - y) * pnorm(z) + s * dnorm(z)
}
sigma_eps <- self$nug * self$s2_hat
sigma_eps2 <- sigma_eps^2
augterm <- function(s2) {
  1 - sigma_eps / sqrt(s2 + sigma_eps2)
}

oout <- optim(
  # par=self$X[star_star_index,],
  par=self$X[star_star_index,] + rnorm(4),
  # par=colMeans(self$X),
  fn=function(x) {
    print(x)
    pred <- self$predict(x, se.fit = T)
    val <- - EI(pred$mean, pred$se) * augterm(pred$s2)
    print(val)
    val
  },
  lower=apply(self$X, 2, min),
  upper=apply(self$X, 2, max),
  method=if (self$D==1) {'Brent'} else {"L-BFGS-B"}
)
oout
list(
  par=oout$par,
  value=oout$value
)



maxAugEI <- function(self) {

  # self <- gp1
  pred_X <- self$predict(self$X, se.fit = T)
  u_X <- -pred_X$mean - pred_X$se
  star_star_index <- which.max(u_X)
  f <- pred_X$mean[star_star_index]
  EI <- function(y, s) {
    z <- (f - y) / s
    (f - y) * pnorm(z) + s * dnorm(z)
  }
  sigma_eps <- self$nug * self$s2_hat
  sigma_eps2 <- sigma_eps^2
  augterm <- function(s2) {
    1 - sigma_eps / sqrt(s2 + sigma_eps2)
  }

  oout <- optim(
    # par=self$X[star_star_index,],
    par=self$X[star_star_index,] + rnorm(4),
    # par=colMeans(self$X),
    fn=function(x) {
      print(x)
      pred <- self$predict(x, se.fit = T)
      val <- - EI(pred$mean, pred$se) * augterm(pred$s2)
      print(val)
      val
    },
    lower=apply(self$X, 2, min),
    upper=apply(self$X, 2, max),
    method=if (self$D==1) {'Brent'} else {"L-BFGS-B"}
  )
  oout
  list(
    par=oout$par,
    value=oout$value
  )
}
