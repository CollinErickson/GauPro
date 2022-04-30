n <- 100
X <- matrix(runif(n), ncol=1)
y <- c(3.2*X-1.4) + rnorm(n)
gpld <- GauPro_kernel_model$new(X, y, kernel='m52')
gpld$plot1D()
gpld$EI(1)
curve(gpld$EI(matrix(x, ncol=1)))
gpld$maxqEI(5, 'pred')

n <- 10
X <- matrix(runif(n), ncol=1)
# y <- c(-3.2*(X-.5)^2-1.4) + rnorm(n,0,.01)
f <- function(X) {c(-3.2*(X-.5)^2-1.4) + rnorm(length(X),0,1e-2)}
y <- f(X)
gpld <- GauPro_kernel_model$new(X, y, kernel='m32')
gpld$plot1D()
gpld$EI(1)
curve(gpld$EI(matrix(x, ncol=1)))
gpld$maxqEI(5, 'pred')
gpld$maxqEI(5, 'CL')
xEI <- gpld$maxqEI(5, 'pred')
gpld$update(Xnew=xEI, Znew=f(xEI))
gpld$plot1D()
plot(gpld$X)


# Knowledge gradient

n <- 10
X <- matrix(runif(n), ncol=1)
# y <- c(-3.2*(X-.5)^2-1.4) + rnorm(n,0,.01)
f <- function(X) {c(-3.2*(X-.5)^2-1.4) + rnorm(length(X),0,1e-2)}
y <- f(X)
gpkg <- GauPro_kernel_model$new(X, y, kernel='m32')
gpkg$plot1D()
# Find current max
gpkgmax <- optim(par=gpkg$X[which.max(gpkg$Z)[1],],
                 fn=function(xx) {-gpkg$pred(xx)}, method='Brent', lower=0, upper=1)
gpkgmax
# Calculate knowledge gradient at xkg
xkg <- .6
# Sample at xkg
xkgpred <- gpkg$pred(xkg, se.fit = T)
xkgpred
nsamps <- 5
xkgsamps <- qnorm(((1:nsamps)-.5)/nsamps, xkgpred$mean, xkgpred$se)
kgs <- rep(NA, nsamps)
for (i in 1:nsamps) {
  xkgsamp <- xkgsamps[i]
  # xkgsamp <- rnorm(1, xkgpred$mean, xkgpred$se)
  # Add samp to mod
  gpkgclone <- gpkg$clone(deep=TRUE)
  gpkgclone$update(Xnew=xkg, Znew=xkgsamp, no_update = TRUE)
  gpkgclone$plot1D()
  # Find clone max after adding sample
  gpkgmaxclone <- optim(par=gpkgclone$X[which.max(gpkgclone$Z)[1],],
                        fn=function(xx) {-gpkgclone$pred(xx)}, method='Brent', lower=0, upper=1)
  gpkgmaxclone
  gpkgmaxclone$value - gpkgmax$value
  kgs[i] <- gpkgmaxclone$value - gpkgmax$value
}
kgs

gpkg <- GauPro_kernel_model$new(X, y, kernel='m32')
gpkg$KG(.6)
curve(sapply(x, function(x) gpkg$KG(x)), n=11)
system.time(curve(sapply(x, function(x) gpkg$KG(x)), n=11))
