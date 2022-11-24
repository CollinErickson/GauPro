# 2 dicrete inputs, 2 cts inputs

f <- function(a, b, c, d) {
  -1e-3*a^2*b^2*a + a*ifelse(c==1,1,0) + ifelse(d==1,1,2) +
    rnorm(length(a),0,1e0)
}
n <- 133
library(dplyr)
Xdf <- bind_cols(
  a=runif(n,6,8),
  b=runif(n,-8,-2),
  c=sample(1:2,n,T),
  d=sample(1:4,n,T)
)
Xmat <- as.matrix(Xdf)
Xmat
# y <- apply(Xmat, 1, f)
y <- f(Xmat[,1],Xmat[,2],Xmat[,3],Xmat[,4])

system.time({
  gp1 <- GauPro_kernel_model$new(
    X=Xmat, Z=y,
    kernel=IgnoreIndsKernel$new(ignoreinds = 3:4, Gaussian$new(D=2))
  )
})
gp1
summary(gp1)
gp1$pred(c(7,-4, 1, 1), se.fit = T)
gp1$nug
gp1$kernel$k(Xmat)

# With 0 restarts
# 23.8,27.9,43.6,19.4 sec before Rcpp on latent
# 8.1, 6.7,10.1 after adding Rcpp
# With 5 restarts
# 146.7,99.6 sec without Rcpp
# 126,93 sec with Rcpp
system.time({
  gp2 <- GauPro_kernel_model$new(
    X=Xmat, Z=y,
    restarts=0, track_optim = T, verbose=5,
    kernel=IgnoreIndsKernel$new(ignoreinds = 3:4, Gaussian$new(D=2)) *
      LatentFactorKernel$new(D=4, nlevels = 2, latentdim = 1, xindex = 3) *
      LatentFactorKernel$new(D=4, nlevels = 4, latentdim = 2, xindex = 4)
  )
})
summary(gp2)
gp2$plotmarginal()
gp2$plotmarginalrandom()
gp2$plot_track_optim()
gp2
gp2$pred(c(7,-4, 1, 1), se.fit = T)
gp2$pred(c(7,-4, 2, 1), se.fit = T)
gp2$pred(c(7,-4, 1, 1), se.fit = T)
gp2$pred(c(7,-4, 1, 2), se.fit = T)
gp2$pred(c(7,-4, 1, 3), se.fit = T)
gp2$pred(c(7,-4, 1, 4), se.fit = T)
gp2$plotLOO()
gp2$nug
gp2$kernel$k(Xmat)

debugonce(gp2$kernel$k2$k)
gp2$kernel$k(Xmat)

# With 0 restarts
# 23.8,27.9,43.6,19.4 sec before Rcpp on latent
# 8.1, 6.7,10.1 after adding Rcpp
# With 5 restarts
# 146.7,99.6 sec without Rcpp
# 126,93 sec with Rcpp
gp2m2 <- profvis::profvis(interval=.01, {
  gp2 <- GauPro_kernel_model$new(
    X=Xmat, Z=y,
    restarts=0,
    kernel=IgnoreIndsKernel$new(ignoreinds = 3:4, Gaussian$new(D=2)) *
      LatentFactorKernel$new(D=4, nlevels = 2, latentdim = 1, xindex = 3) *
      LatentFactorKernel$new(D=4, nlevels = 4, latentdim = 2, xindex = 4)
  )
})

# Use formula
gp2f <- GauPro_kernel_model$new(data=Xdf, y ~ a+b+factor(c)+factor(d))
gp2f
gp2f$kernel
gp2f$plot()
gp2f$summary()

# Only factors
system.time({
  gp3 <- GauPro_kernel_model$new(
    X=Xmat[,3:4], Z=y,
    restarts=0, track_optim = T, verbose=5,
    kernel=LatentFactorKernel$new(D=2, nlevels = 2, latentdim = 1, xindex = 1) *
      LatentFactorKernel$new(D=2, nlevels = 4, latentdim = 2, xindex = 2)
  )
})
gp3$maxEI()








# gp4: 2 disc, 2 cts, 2 int ----
f <- function(a, b, c, d, e, f) {
  -1e-3*a^2*b^2*a + a*ifelse(c==1,1,0) + ifelse(d==1,1,2) + rnorm(length(a),0,1e-1) +
    a*e/1e3 + .3*f + sin(f)*e/1e3
}
n <- 133
library(dplyr)
Xdf <- bind_cols(
  a=runif(n,6,8),
  b=runif(n,-8,-2),
  c=sample(1:2,n,T),
  d=sample(1:4,n,T),
  e=sample(1:10000, n, T),
  f=sample(c(1,3,5,7,9),n,T)
)
Xmat <- as.matrix(Xdf)
Xmat
# y <- apply(Xmat, 1, f)
y <- f(Xmat[,1],Xmat[,2],Xmat[,3],Xmat[,4],Xmat[,5],Xmat[,6])

system.time({
  gp4 <- GauPro_kernel_model$new(
    X=Xmat, Z=y,
    kernel=IgnoreIndsKernel$new(ignoreinds = 3:4, Gaussian$new(D=4)) *
      LatentFactorKernel$new(D=4, nlevels = 2, latentdim = 1, xindex = 3) *
      LatentFactorKernel$new(D=4, nlevels = 4, latentdim = 2, xindex = 4),
    nug.max=.1
  )
})
gp4
plot(gp4$pred(Xmat), y); abline(a=0,b=1, col=2)
gp4$plotLOO()
gp4$maxEI()
gp4$maxEI(discreteinputs = list('5'=1:1e4, '6'=c(1,3,5,7,9)))
gp4$maxEIwithfactorsordiscrete2(discreteinputs = list('5'=1:1e4, '6'=c(1,3,5,7,9)))

