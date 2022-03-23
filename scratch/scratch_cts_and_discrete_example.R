# 2 dicrete inputs, 2 cts inputs

f <- function(a, b, c, d) {
  -1e-3*a^2*b^2*a + ifelse(d==1,1,2) + rnorm(length(a),0,1e-1)
}
n <- 33
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

gp1 <- GauPro_kernel_model$new(
  X=Xmat, Z=y,
  kernel=IgnoreIndsKernel$new(ignoreinds = 3:4, Gaussian$new(D=2))
)
gp1
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
    restarts=0,
    kernel=IgnoreIndsKernel$new(ignoreinds = 3:4, Gaussian$new(D=2)) *
      LatentFactorKernel$new(D=4, nlevels = 2, latentdim = 1, xindex = 3) *
      LatentFactorKernel$new(D=4, nlevels = 4, latentdim = 2, xindex = 4)
  )
})
gp2
gp2$pred(c(7,-4, 1, 1), se.fit = T)
gp2$nug
gp2$kernel$k(Xmat)

debugonce(gp2$kernel$k2$k)
gp2$kernel$k(Xmat)
