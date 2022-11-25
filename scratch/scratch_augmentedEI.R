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
  # Get preds at existing points, calculate best
  pred_X <- self$predict(self$X, se.fit = T)
  u_X <- -pred_X$mean - pred_X$se
  star_star_index <- which.max(u_X)

  f <- pred_X$mean[star_star_index]

  # Func to calculate EI
  EI <- function(y, s) {
    z <- (f - y) / s
    (f - y) * pnorm(z) + s * dnorm(z)
  }
  # Calculate "augmented" term
  sigma_eps <- self$nug * self$s2_hat
  sigma_eps2 <- sigma_eps^2
  augterm <- function(s2) {
    1 - sigma_eps / sqrt(s2 + sigma_eps2)
  }

  if (F) {
    # derivative:
    # dAugEI_dx = ei * daugterm_dx + dei_dx * augterm
    # daug_dx = .5*seps / (s2 + seps2)^1.5 * ds2_dx
    # dEI_dx = -dy_dx*pnorm(z) + (f-y)*dnorm(z)*dz_dx + ds_dx*dnorm(z) + s*ddnormz_dz(z)*dz_dx
    # ddnormz_dz = dnorm(z) * (-2)*z
    # dz_dx = -dy_dx/s + (f-y)(-1/s^2)*ds_dx
    # ds_dx = .5/s * ds2_dx

    x <- .8
    predx <- self$pred(x, se=T)
    y <- predx$mean
    s <- predx$se
    s2 <- predx$s2
    ds2_dx <- self$gradpredvar(x) # GOOD
    ds_dx <- .5/s * ds2_dx # GOOD
    z <- (f - y) / s
    dy_dx <- self$grad(x) # GOOD
    dz_dx <- -dy_dx / s + (f - y) * (-1/s2) * ds_dx # GOOD
    ddnormz_dz <- -dnorm(z) * z # GOOD
    daug_dx = .5*sigma_eps / (s2 + sigma_eps2)^1.5 * ds2_dx # GOOD
    dEI_dx = -dy_dx*pnorm(z) + (f-y)*dnorm(z)*dz_dx + ds_dx*dnorm(z) + s*ddnormz_dz*dz_dx #GOOD
    numDeriv::grad(function(x) {pr <- self$pred(x,se=T);( EI(pr$mean,pr$se))}, x)
    dAugEI_dx = EI(y, s) * daug_dx + dEI_dx * augterm(s2)
    numDeriv::grad(function(x) {pr <- self$pred(x,se=T);( EI(pr$mean,pr$se)*augterm(pr$s2))}, x)
  }

  # Optimize
  oout <- optim(
    # par=self$X[star_star_index,],
    par=self$X[star_star_index*1+0,], #+ rnorm(ncol(self$X)),
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
    value=-oout$value
  )
}
maxAugEI(gp2)
maxAugEI(gp)


AugEI <- function(self, x, minimize=FALSE, eps=0,
                  return_grad=F, f=NULL) {
  stopifnot(length(minimize)==1, is.logical(minimize))
  stopifnot(length(eps)==1, is.numeric(eps), eps >= 0)
  if (is.matrix(x)) {
    stopifnot(ncol(x) == ncol(self$X))
  } else if (is.vector(x)) {
    stopifnot(length(x) == ncol(self$X))
  } else if (is.data.frame(x) && !is.null(self$formula)) {
    # Fine here, will get converted in predict
  } else {
    stop(paste0("bad x in EI, class is: ", class(x)))
  }

  if (is.null(f)) {
    # Get preds at existing points, calculate best
    pred_X <- self$predict(self$X, se.fit = T)
    if (minimize) {
      u_X <- -pred_X$mean - pred_X$se
      star_star_index <- which.max(u_X)
    } else {
      warning("AugEI must minimize for now")
      u_X <- +pred_X$mean + pred_X$se
      star_star_index <- which.max(u_X)
    }

    f <- pred_X$mean[star_star_index]

  } else {
    stopifnot(is.numeric(f), length(f) == 1)
  }


  predx <- self$pred(x, se=T)
  y <- predx$mean
  s <- predx$se
  s2 <- predx$s2

  z <- (f - y) / s
  EI <- (f - y) * pnorm(z) + s * dnorm(z)


  # Calculate "augmented" term
  sigma_eps <- self$nug * self$s2_hat
  sigma_eps2 <- sigma_eps^2
  Aug <- 1 - sigma_eps / sqrt(s2 + sigma_eps2)
  AugEI <- Aug * EI

  if (return_grad) {
    # x <- .8
    ds2_dx <- self$gradpredvar(x) # GOOD
    ds_dx <- .5/s * ds2_dx # GOOD
    # z <- (f - y) / s
    dy_dx <- self$grad(x) # GOOD
    dz_dx <- -dy_dx / s + (f - y) * (-1/s2) * ds_dx # GOOD
    ddnormz_dz <- -dnorm(z) * z # GOOD
    daug_dx = .5*sigma_eps / (s2 + sigma_eps2)^1.5 * ds2_dx # GOOD
    dEI_dx = -dy_dx*pnorm(z) + (f-y)*dnorm(z)*dz_dx + ds_dx*dnorm(z) + s*ddnormz_dz*dz_dx #GOOD
    # numDeriv::grad(function(x) {pr <- self$pred(x,se=T);( EI(pr$mean,pr$se))}, x)
    dAugEI_dx = EI * daug_dx + dEI_dx * Aug
    # numDeriv::grad(function(x) {pr <- self$pred(x,se=T);( EI(pr$mean,pr$se)*augterm(pr$s2))}, x)
    return(list(
      AugEI=AugEI,
      grad=dAugEI_dx
    ))
  }
  AugEI
}
# 1D
u <- .8
AugEI(self, u)
AugEI(self, u, return_grad = T)
numDeriv::grad(function(u)AugEI(self,u), u)

u <- matrix(runif(3), ncol=1)
AugEI(self, u)
AugEI(self, u, return_grad = T)
numDeriv::grad(function(u)AugEI(self,u), .8)
curve(AugEI(self, x))
curve(self$EI(matrix(x, ncol=1),minimize = T))

# 2D

u <- matrix(c(.5,0), ncol=2)
AugEI(self, u)
AugEI(self, u, return_grad = T)
numDeriv::grad(function(u)AugEI(self,u), u)

v <- rbind(matrix(runif(6), ncol=2), u)
AugEI(self, v)
AugEI(self, v, return_grad = T)
numDeriv::grad(function(u)AugEI(self,u), .8)
curve(AugEI(self, x))
curve(self$EI(matrix(x, ncol=1),minimize = T))

maxAugEI <- function() {

}
