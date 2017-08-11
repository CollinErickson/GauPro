# Check grad
set.seed(0)
n <- 20
x <- matrix(seq(0,1,length.out = n), ncol=1)
f <- Vectorize(function(x) {sin(2*pi*x) + .5*sin(4*pi*x) +rnorm(1,0,.3) + 10*x})
y <- 123 + f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
kernels <- list(Gaussian$new(1), Exponential$new(1), Matern32$new(0), Matern52$new(0), Periodic$new(.1, .2, alpha_lower=.1, p_lower=.1), RatQuad$new(.1, alpha=2.4), Exponential$new(1)+ Matern32$new(0), Exponential$new(1, beta_lower=-1)* Matern32$new(0, beta_lower=-1) )
trends <- list(trend_LM$new(D=1), trend_0$new())
kernel <- kernels[[8]]
trend <- trends[[2]]
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=kernel, trend=trend, parallel=FALSE, verbose=10, nug.est=T, restarts=1)
gp$cool1Dplot()

XX <- matrix(c(.5, .6), ncol=1)
numDeriv::grad(gp$predict, XX)
gp$grad(XX)


# check 2D

set.seed(0)
n <- 100
x <- matrix(runif(2*40), ncol=2)
f <- function(x) {sin(2*pi*x[1]) + .5*sin(4*pi*x[2]) + x[2]^.7 + rnorm(1,0,1)}
y <- apply(x, 1, f) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
kernels <- list(Gaussian$new(0:1), Exponential$new(0:1), Matern32$new(0:1), Matern52$new(0:1), Periodic$new(alpha=.1, p=c(.2,.3), alpha_lower=.1, p_lower=c(.1,.1), p_upper=c(10,10)), RatQuad$new(c(.1,.2), alpha=2.), Exponential$new(0:1)+Matern32$new(0:1), Exponential$new(0:1)*Matern32$new(0:1))
trends <- list(trend_LM$new(D=2))
kernel <- kernels[[8]]
trend <- trends[[1]]
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=kernel, nug=.1, trend=trend, parallel=FALSE, verbose=10, nug.est=T, restarts=0)
# ContourFunctions::cf(gp$predict, pts=x, batchmax=300) #cool1Dplot()

XX <- matrix(c(.5, .6, .7, .8, .9, .95), ncol=2)
# numDeriv::grad(gp$predict, XX)
gp$grad(XX)
t(apply(XX, 1, function(xx) numDeriv::grad(gp$predict, xx)))
