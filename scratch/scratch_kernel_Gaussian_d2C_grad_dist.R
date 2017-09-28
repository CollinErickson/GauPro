# Check dC_dx
# For Gaussian_beta kernel
set.seed(0)
n <- 20
x <- matrix(seq(0,1,length.out = n), ncol=1)
f <- Vectorize(function(x) {sin(2*pi*x) + .5*sin(4*pi*x) +rnorm(1,0,.3)}+10*x)
y <- 123 + f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian$new(1), parallel=FALSE, verbose=10, nug.est=T, s2_est=F)
gp$cool1Dplot()

# Following two should be equal
gp$kernel$dC_dx(XX=matrix(.343), X=gp$X)[1,1,]
sapply(1:20,
       function(i) numDeriv::grad(function(xx) {gp$kernel$k(matrix(xx), x[i,])}, x=.343)
)

# This should work too
gp$kernel$dC_dx(XX=matrix(.343), X=matrix(.56))[1,1,]
numDeriv::grad(function(xx) {gp$kernel$k(matrix(xx), matrix(.56))}, x=.343)

# Check 2D
set.seed(0)
n <- 30
x <- lhs::maximinLHS(n=n, k=2)
f <- function(x) {sin(2*pi*x[1]) + .5*sin(4*pi*x[1]) +rnorm(1,0,.03) + x[2]^2}
y <- apply(x, 1, f) #f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian, parallel=FALSE, verbose=10, nug.est=T)
ContourFunctions::cf(gp$predict, pts=x)
gp$kernel$dC_dx(XX=matrix(c(.1,.343),ncol=2), X=matrix(c(.34,.56),ncol=2))
numDeriv::grad(function(xx) {gp$kernel$k(matrix(xx,ncol=2), matrix(c(.34,.56),ncol=2))}, x=c(.1,.343))



# Check d2C_dx2
# For Gaussian_beta kernel
set.seed(0)
n <- 20
x <- matrix(seq(0,1,length.out = n), ncol=1)
f <- Vectorize(function(x) {sin(2*pi*x) + .5*sin(4*pi*x) +rnorm(1,0,.3)}+10*x)
y <- 123 + f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian$new(1), parallel=FALSE, verbose=10, nug.est=T, s2_est=F)
gp$cool1Dplot()

# Following two should be equal
gp$kernel$d2C_dx2(XX=matrix(.343), X=gp$X)[,1,1,]
sapply(1:20,
       function(i) numDeriv::hessian(function(xx) {gp$kernel$k(matrix(xx), x[i,])}, x=.343)
)

# This should work too
gp$kernel$d2C_dx2(XX=matrix(.343), X=matrix(.56))[,1,1,]
numDeriv::hessian(function(xx) {gp$kernel$k(matrix(xx), matrix(.56))}, x=.343)

# Check 2D
set.seed(0)
n <- 30
x <- lhs::maximinLHS(n=n, k=2)
f <- function(x) {sin(2*pi*x[1]) + .5*sin(4*pi*x[1]) +rnorm(1,0,.03) + x[2]^2}
y <- apply(x, 1, f) #f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian, parallel=FALSE, verbose=10, nug.est=T)
ContourFunctions::cf(gp$predict, pts=x)
gp$kernel$d2C_dx2(XX=matrix(c(.1,.343),ncol=2), X=matrix(c(.34,.56),ncol=2))
numDeriv::hessian(function(xx) {gp$kernel$k(matrix(xx,ncol=2), matrix(c(.34,.56),ncol=2))}, x=c(.1,.343))

# Want cross d2C_dudv
gp$kernel$d2C_dudv(XX=matrix(c(.1,.343),ncol=2), X=matrix(c(.34,.56),ncol=2))
# tt <- numDeriv::genD(function(xx) {gp$kernel$k(matrix(xx,ncol=2), matrix(c(.34,.56),ncol=2))}, x=c(.1,.343))
numDeriv::jacobian(
  function(x2) {
    numDeriv::grad(
      function(xx) {
        gp$kernel$k(matrix(xx,ncol=2), matrix(x2,ncol=2))
      },
      #x=matrix(x2,ncol=2))
      x=matrix(c(.34,.56),ncol=2))
}, x=c(.1,.343))


gp$grad_dist(XX=matrix(c(.1,.343),ncol=2))
eps <- 1e-5
ts <- gp$sample(XX = matrix(c(.1+eps, .343+eps, .1+eps, .343-eps, .1-eps, .343+eps, .1-eps, .343-eps), ncol=2, byrow = T), n = 1e6)
ts2 <- apply(ts, 1, function(a) {1/eps/2 * c(a[2]-a[4], a[3]-a[4])})
rbind(rowMeans(ts2), apply(ts2, 1, sd))
rbind(rowMeans(ts2), cov(t(ts2)))


# Check g2 norm mean
xx <- c(.22,.554)
xx <- c(.1,.343)
tg <- gp$grad_dist(XX=matrix(xx,ncol=2))
m <- tg$mean[1,]
Sigma <- tg$cov[1,,]
SigmaInv <- solve(Sigma)
SigmaInvRoot <- expm::sqrtm(SigmaInv)
th <- tg$cov[1,,]
eth <- eigen(th)
P <- t(eth$vectors)
lambda <- eth$values
testthat::expect_equal(t(P) %*% diag(eth$values) %*% (P), Sigma) # Should be equal
b <- P %*% SigmaInvRoot %*% m
c(sum(b^2 * lambda)+length(m), 4*sum(b^2 * lambda^2)+2*length(m))
eps <- 1e-4
gsampraw <- gp$sample(XX = matrix(c(xx[1]+eps, xx[2]+eps, xx[1]+eps, xx[2]-eps, xx[1]-eps, xx[2]+eps, xx[1]-eps, xx[2]-eps), ncol=2, byrow = T), n = 1e6)
gsamp <- apply(gsampraw, 1, function(a) {1/eps/2 * c(a[2]-a[4], a[3]-a[4])})
g2samp <- apply(gsampraw, 1, function(a) {sum((1/eps/2 * c(a[2]-a[4], a[3]-a[4]))^2)})
tg
rbind(rowMeans(gsamp), cov(t(gsamp)))
c(mean(g2samp), var(g2samp))
X <- 4:5
Y <- SigmaInvRoot %*% X
Z <- Y - SigmaInvRoot %*% m
U <- P %*% Z
testthat::expect_equal(sum(X^2), sum(lambda * (U + b) ^ 2))

# Not working, so check values
Xsamp <- gsamp
Ysamp <- apply(Xsamp, 2, function(a) {SigmaInvRoot %*% a})
Zsamp <- sweep(Ysamp,1, SigmaInvRoot %*% m)
Usamp <- apply(Zsamp, 2, function(zz) {P %*% zz})
apply(Usamp, 1, summary)
cov(t(Usamp))

# Make sure it scales linearly
xr <- 2^ (3:9)
yr <- sapply(xr, function(xi){print(xi);system.time(gp$grad_norm2_dist(matrix(c(.2,.3),ncol=2, nrow=xi)))[3]})
plot(xr, yr)

# Check plot
set.seed(30)
n <- 30
x <- lhs::maximinLHS(n=n, k=2)
f <- TestFunctions::banana #function(x) {sin(2*pi*x[1]) + .5*sin(4*pi*x[1]) +rnorm(1,0,.03) + x[2]^2}
y <- apply(x, 1, f) #f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian, parallel=FALSE, verbose=10, nug.est=T)
ContourFunctions::cf(gp$predict, batchmax=Inf, pts=gp$X)
ContourFunctions::cf(function(a)gp$grad_norm2_dist(a)$mean, batchmax=Inf, pts=gp$X, n=100)


# See why sqrtm gives message, see if other way
Kx <- gp$kernel$k(x)
expm::sqrtm(Kx)
e <- eigen(Kx)
V <- e$vectors
summary(c(V %*% diag(e$values) %*% t(V) - Kx))
B <- V %*% diag(sqrt(e$values)) %*% t(V)
mysqrt = function(mat, symmetric) {
  e <- eigen(mat, symmetric=symmetric)
  V <- e$vectors
  B <- V %*% diag(sqrt(e$values)) %*% t(V)
  B
}
summary(c(mysqrt(Kx) - expm::sqrtm(Kx)))
summary(c(mysqrt(Kx, T) - expm::sqrtm(Kx)))
m1 <- mysqrt(Kx)
m2 <- mysqrt(Kx,T)
m3 <- expm::sqrtm(Kx)
mean((c(m1 - Kx)^2))
mean((c(m2 - Kx)^2))
mean((c(m3 - Kx)^2))
microbenchmark::microbenchmark(mysqrt(Kx), mysqrt(Kx, T), expm::sqrtm(Kx))
