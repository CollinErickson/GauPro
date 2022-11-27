# EI with high noise. Gives too much improvement.
d <- 1
n <- 20
# x <- runif(n)
x <- seq(0,1,l=n)
x <- c(runif(n/2,0,.3), runif(5,.3,.7), runif(n/2-5,.7,1))
y <- sin(2*pi*x) + rnorm(length(x), 0, .2)
plot(x, y)
# gp1 <- GauPro_kernel_model$new(y ~ x)
gp1 <- GauPro_kernel_model$new(x, y, kernel='m52')
gp1$cool1Dplot()
gp1$plot1D()
curve(gp1$pred(x, se.fit = T)$se, ylim=c(0,.3))
curve(gp1$pred(x, se.fit = T, mean_dist = T)$se, add=T, col=2)
# EI is way too big
gp1$maxEI(minimize = T)
curve(sapply(x, function(x)gp1$EI(x, minimize = T)))


# What should it be?
# If you add the point.
# Sample over dist.
# Add that sample to the model.
# Find new mean prediction at that point.
# KG runs optimization after that.


EI2 <- function(self, x, minimize=FALSE) {
  self <- gp1
  predX <- self$pred(self$X)

  x <- .35
  X <- self$X
  Z <- self$Z
  u <- x
  a <- x

  mu_a <- self$trend$Z(a)
  mu_X <- c(self$trend$Z(X))
  KX <- self$K
  KXinv <- self$Kinv
  KuX <- matrix(self$kernel$k(u, X), nrow=1)
  KXu <- t(KuX)
  KaX <- self$kernel$k(a, X)
  KXa <- t(KaX)
  Kau <- self$kernel$k(a, u)
  Ka <- self$kernel$k(a)
  Ku <- self$kernel$k(u) + self$s2_hat*self$nug
  Ma <- c(KaX, Kau)
  Mcinv <- rbind(cbind(KX, KXu),
                 cbind(KuX, Ku))
  Mc <- solve(Mcinv)
  KugivenX <- Ku - KuX %*% KXinv %*% KXu

  self$pred(a, se=T, mean_dist = T)

  mu_m <- c(mu_a + Ma %*% Mc %*% rbind(Z-mu_X, c(KuX%*%KXinv%*%(Z-mu_X))))
  mu_m

  cov_m <- matrix(Ma, nrow=1) %*% Mc %*%
    rbind(cbind(matrix(0, ncol=ncol(KX)+1, nrow=nrow(KX))),
          cbind(matrix(0, ncol=ncol(KX),   nrow=1),
                KugivenX)) %*%
    Mc %*% matrix(Ma, ncol=1)
  cov_m
  # This is right! Matches simulation below.


  # Sample to estimate dist
  ns <- 1e3
  ms <- rep(0, ns)
  for (i in 1:ns) {
    gpi <- self$clone(deep=T)
    zi <- gpi$sample(u)
    gpi$update(Xnew=u, Znew=zi, no_update = T)
    ms[i] <- gpi$pred(u)
  }
  summary(ms)
  var(ms)
  hist(ms, freq=F)
  curve(dnorm(x, mu_m, sqrt(c(cov_m))), add=T, col=2)
}




CorrectedEI <- function(self, x, minimize=FALSE, eps=0,
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
    pred_X <- self$predict(self$X, se.fit = F)
    if (minimize) {
      # u_X <- -pred_X$mean - pred_X$se
      star_star_index <- which.min(pred_X)
    } else {
      # warning("AugEI must minimize for now")
      # u_X <- +pred_X$mean + pred_X$se
      star_star_index <- which.max(pred_X)
    }

    f <- pred_X[star_star_index]

  }
  stopifnot(is.numeric(f), length(f) == 1)



  # predx <- self$pred(x, se=T)
  # y <- predx$mean
  # s <- predx$se
  # s2 <- predx$s2

  minmult <- if (minimize) {1} else {-1}

  x <- matrix(seq(0,1,l=131), ncol=1)
  u <- x
  X <- self$X
  mu_u <- self$trend$Z(u)
  Ku.X <- self$kernel$k(u, X)
  mu_X <- self$trend$Z(X)
  Ka <- self$kernel$k(u)
  Ku <- Ka + self$nug * self$s2_hat
  Ku_given_X <- Ku - Ku.X %*% self$Kinv %*% t(Ku.X)

  y <- mu_u + Ku.X %*% self$Kinv %*% (self$Z - mu_X)
  s2 <- (Ku_given_X - self$nug*self$s2_hat) ^ 2 / (Ku_given_X)
  if (ncol(s2) > 1.5) {s2 <- diag(s2)}
  s <- sqrt(s2)

  # int from f to Inf: (x-f) p(x) dx


  z <- (f - y) / s * minmult
  CorEI <- (f - y) * minmult * pnorm(z) + s * dnorm(z)
  tdf <- 3
  CorEIt <- (f - y) * minmult * pt(z,tdf) + s * dt(z,tdf)
  plot(x, CorEI)
  plot(x, s, ylim=c(0,.3))
  points(x, self$pred(x, se=T)$se,col=2)
  points(x, self$pred(x, se=T, mean_dist = T)$se,col=3)
  cbind(x, y, s, z, CorEI=CorEI, EIt=(f - y) * minmult * pt(z,3) + s * dt(z, 3))


  # # Calculate "augmented" term
  # sigma_eps <- self$nug * self$s2_hat
  # sigma_eps2 <- sigma_eps^2
  # Aug <- 1 - sigma_eps / sqrt(s2 + sigma_eps2)
  # AugEI <- Aug * EI

  if (F && return_grad) {
    # x <- .8
    ds2_dx <- self$gradpredvar(x) # GOOD
    ds_dx <- .5/s * ds2_dx # GOOD
    # z <- (f - y) / s
    dy_dx <- self$grad(x) # GOOD
    dz_dx <- -dy_dx / s + (f - y) * (-1/s2) * ds_dx # GOOD
    dz_dx <- dz_dx * minmult
    ddnormz_dz <- -dnorm(z) * z # GOOD
    daug_dx = .5*sigma_eps / (s2 + sigma_eps2)^1.5 * ds2_dx # GOOD
    dEI_dx = minmult * (-dy_dx*pnorm(z) + (f-y)*dnorm(z)*dz_dx) +
      ds_dx*dnorm(z) + s*ddnormz_dz*dz_dx #GOOD
    # numDeriv::grad(function(x) {pr <- self$pred(x,se=T);( EI(pr$mean,pr$se))}, x)
    dAugEI_dx = EI * daug_dx + dEI_dx * Aug
    # numDeriv::grad(function(x) {pr <- self$pred(x,se=T);( EI(pr$mean,pr$se)*augterm(pr$s2))}, x)
    return(list(
      AugEI=AugEI,
      grad=dAugEI_dx
    ))
  }
  CorEI
}
