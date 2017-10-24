
n <- 12
x <- matrix(seq(0,1,length.out = n), ncol=1)
f <- function(x) {sin(2*pi*x)}
y <- f(x) + rnorm(n,0,1e-2)

gp <- GauPro_Gauss_LOO$new(X=x, Z=y, parallel=FALSE)

curve(f, ylim=c(-2,2))
curve(gp$predict(x), add=T, col=2)
curve(gp$predict(x) + 2 * gp$predict(x,se=T)$se, add=T, col=3)
curve(gp$predict(x) - 2 * gp$predict(x,se=T)$se, add=T, col=3)
points(x, y)

curve(gp$tmod$predict(x))
points(gp$X, gp$tmod$Z)

# Try function with big errors at middle
n <- 13
x <- runif(n)
f <- function(x) {1 / sqrt(1 + exp(-25*(x-.5)))}
y <- f(x) + rnorm(n,0,1e-3)

gp <- GauPro_Gauss_LOO$new(X=x, Z=y, parallel=FALSE)

curve(f, ylim=c(-.05,1.1))
curve(gp$predict(x), add=T, col=2)
curve(gp$predict(x) + 2 * gp$predict(x,se=T)$se, add=T, col=3)
curve(gp$predict(x) - 2 * gp$predict(x,se=T)$se, add=T, col=3)
points(x, y)


# Try function with big errors at middle
n <- 20
x <- runif(n)
f <- function(x) {(2*x) %% 1}
y <- f(x) + rnorm(n,0,1e-1)

gp <- GauPro_Gauss_LOO$new(X=x, Z=y, parallel=FALSE)

curve(f, ylim=c(-.05,1.1))
curve(gp$predict(x), add=T, col=2)
curve(gp$predict(x) + 2 * gp$predict(x,se=T)$se, add=T, col=3)
curve(gp$predict(x) - 2 * gp$predict(x,se=T)$se, add=T, col=3)
points(x, y)
curve(gp$tmod$predict(x))


# Try function that has bigger errors to right
n <- 20
x <- runif(n)
f <- function(x) {x}
y <- f(x) + rnorm(n,0,x/1)

gp <- GauPro_Gauss_LOO$new(X=x, Z=y, parallel=FALSE)

curve(f, ylim=c(-.05,1.1))
curve(gp$predict(x), add=T, col=2)
curve(gp$predict(x) + 2 * gp$predict(x,se=T)$se, add=T, col=3)
curve(gp$predict(x) - 2 * gp$predict(x,se=T)$se, add=T, col=3)
points(x, y)
curve(gp$tmod$predict(x))
points(x, gp$tmod$Z)
cbind(x, y, gp$pred_LOO(se=T))

# Check banana function
n <- 50
x <- lhs::maximinLHS(n=n, k=2)
y <- TestFunctions::banana(x)
gp <- GauPro_Gauss_LOO$new(X=x, Z=y, parallel=FALSE)
ContourFunctions::cf(gp$predict, batchmax=Inf, pts=x)
ContourFunctions::cf(function(x) gp$predict(x, se.fit=T)$se, batchmax=Inf, pts=x)
ContourFunctions::cf(gp$tmod$predict, batchmax=Inf, pts=x)

# Turn LOO off to see difference
gp$use_LOO <- FALSE
ContourFunctions::cf(gp$predict, batchmax=Inf, pts=x)
ContourFunctions::cf(function(x) gp$predict(x, se.fit=T)$se, batchmax=Inf, pts=x)


# Check on banana function
n <- 80
d <- 2
f1 <- TestFunctions::banana#function(x) {abs(sin(2*pi*x[1]))}
X1 <- matrix(runif(n*d),n,d)
Z1 <- apply(X1,1,f1)# + rnorm(n, 0, 1e-3)
gp <- GauPro_Gauss_LOO$new(X=x, Z=y, parallel=FALSE)
nn <- 1e3
XX <- matrix(runif(nn*d),nn,d)
ZZ <- apply(XX, 1, f1)
gp$use_LOO <- T
ZZhat <- gp$predict(XX, se=T)
plot(ZZ, ZZhat$me)
abline(a=0, b=1)
points(ZZ, ZZhat$me + 2*ZZhat$se, col=3)
points(ZZ, ZZhat$me - 2*ZZhat$se, col=2)
((ZZhat$me - ZZ) / ZZhat$se) %>% abs %>% summary
ContourFunctions::cf(gp$predict, pts=X1)



# Check if LOO predictions match actual on banana function,  i.e. check shortcut
n <- 80
d <- 2
f1 <- TestFunctions::banana#function(x) {abs(sin(2*pi*x[1]))}
X1 <- matrix(runif(n*d),n,d)
Z1 <- apply(X1,1,f1) * 9.3 # + rnorm(n, 0, 1e-3)
gp <- GauPro_Gauss_LOO$new(X=X1, Z=Z1, parallel=FALSE)
ContourFunctions::cf(gp$predict, pts=X1)
nn <- 1e3
XX <- matrix(runif(nn*d),nn,d)
ZZ <- apply(XX, 1, f1) * 9.3
gp$use_LOO <- T
ZZhat <- gp$predict(XX, se=T)

ZLOO <- gp$pred_LOO(se=T)
gp2 <- gp$clone(deep=T)
loo_means <- numeric(n)
loo_ses <- numeric(n)
for (i in 1:n) {
  gpi <- gp$clone(deep=T);
  gpi$update(Xall=X1[-i,],Zall=Z1[-i], no_update = TRUE);
  if (T) { #set mu and s2 back to original values
    # This makes differences ~ 1e-15 instead of 1e-4, not sure if it is recommended though
    gpi$s2_hat <- gp$s2_hat
    gpi$mu_hat <- gp$mu_hat
  }
  gpp <- gpi$predict(X1[i,],se=T)
  loo_means[i] <- gpp$me
  loo_ses[i] <- gpp$se
}
cbind(ZLOO$fit, loo_means)
summary(ZLOO$fit - loo_means)
c(gp$mu_hat, gpi$mu_hat)
c(gp$s2_hat, gpi$s2_hat)


# Check banana again
n <- 40
d <- 2
f1 <- TestFunctions::banana#function(x) {abs(sin(2*pi*x[1]))}
X1 <- MaxPro::MaxProLHD(n=n,p=d, itermax = 10)$Design#t(lhs::maximinLHS(n=d,k=n)) #
# X1 <- matrix(runif(n*d),n,d)
Z1 <- apply(X1,1,f1) * 9.3 # + rnorm(n, 0, 1e-3)
gp <- GauPro_Gauss_LOO$new(X=X1, Z=Z1, parallel=FALSE)
ContourFunctions::cf(gp$predict, pts=X1)
gp_noLOO <- gp$clone(deep = T); gp_noLOO$use_LOO <- FALSE
ContourFunctions::cf(gp_noLOO$predict, pts=X1)
ContourFunctions::cf(gp$tmod$predict, pts=X1)
# Plot se's
ContourFunctions::cf(function(x) gp$pred(x, se=T)$se, pts=X1, batchmax=Inf)
ContourFunctions::cf(function(x) gp_noLOO$pred(x, se=T)$se, pts=X1, batchmax=Inf)
# See predicted t values
nn <- 1e3
XX <- matrix(runif(nn*d),nn,d)
ZZ <- apply(XX, 1, f1) * 9.3
ploo <- gp$pred(XX, se=T)
tloo <- (ploo$mean - ZZ) / ploo$se
summary(tloo)
pnoloo <- gp_noLOO$pred(XX, se=T)
tnoloo <- (pnoloo$mean - ZZ) / pnoloo$se
summary(tnoloo)
stripchart(list(LOO=tloo,no_LOO=tnoloo))
ContourFunctions::cf(function(x){z <- apply(x, 1, f1)*9.3; ploo <- gp$pred(x, se=T);tloo <- (ploo$mean - z) / ploo$se;abs(tloo)}, pts=X1, batchmax=Inf)
ContourFunctions::cf(function(x){z <- apply(x, 1, f1)*9.3; ploo <- gp_noLOO$pred(x, se=T);tloo <- (ploo$mean - z) / ploo$se;abs(tloo)}, pts=X1, batchmax=Inf)
qqnorm(tloo)
qqline(tloo)
qqnorm(tnoloo)
qqline(tnoloo)

