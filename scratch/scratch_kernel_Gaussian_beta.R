# Check numerically that gradient is correct for 1D
# For Gaussian_beta kernel
set.seed(0)
n <- 20
x <- matrix(seq(0,1,length.out = n), ncol=1)
f <- Vectorize(function(x) {sin(2*pi*x) + .5*sin(4*pi*x) +rnorm(1,0,.3)}+10*x)
y <- 123 + f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian_beta$new(1), parallel=FALSE, verbose=10, nug.est=T)
gp$cool1Dplot()
numDeriv::grad(func = function(x)gp$deviance(params=x[2:3], nuglog=x[4], trend_params=x[1]), x=c(10, 2,1, -4))
gp$deviance_grad(params = c(2,1), nug.update=T, nuglog=-4, trend_params=10)
numDeriv::grad(func = function(x)gp$deviance(trend_params=x[1], params=x[2:3], nuglog=x[4]), x=c(gp$trend$m, gp$kernel$beta, gp$kernel$logs2, log(gp$nug,10)))
gp$deviance_grad(params = c(gp$kernel$beta, gp$kernel$logs2), trend_params=gp$trend$m, nug.update=T, nuglog=log(gp$nug,10))

# Check dC_dparams
beta <- .6
s2 <- .3
nug <- 1e-4*10
m1 <- (gp$kernel$k(gp$X, beta=beta+1e-6, s2=s2) - gp$kernel$k(gp$X, beta=beta-1e-6, s2=s2)) / 1e-6/2
C_nonug <- gp$kernel$k(gp$X, beta=beta, s2=s2)
C <- C_nonug + diag(s2*nug, nrow(C_nonug))
m2 <- gp$kernel$dC_dparams(params = c(beta, log(s2,10)), X = gp$X, C = C, C_nonug = C_nonug)[[1]][[1]]
c(m1-m2) %>% summary

# Check if not passing C and C_nonug is okay
gp$kernel$dC_dparams(params = c(beta, log(s2,10)), X = gp$X, C = C, C_nonug = C_nonug)[[1]][[1]]
gp$kernel$dC_dparams(params = c(beta, log(s2,10)), X = gp$X, nug=nug)[[1]][[1]]


# Check C_dC_dparams
params <- c(1.2,.8)
nug <- .001
gp$deviance(params=params, nug=nug)
gp$deviance_grad(params=params, nug=nug, nug.update=T)
gp$deviance_fngr(params=params, nug=nug, nug.update=T)
microbenchmark::microbenchmark(sep={gp$deviance(params=params, nug=nug);gp$deviance_grad(params=params, nug=nug, nug.update=T)}, fngr=gp$deviance_fngr(params=params, nug=nug, nug.update=T))

# Check 2D
set.seed(0)
n <- 30
x <- lhs::maximinLHS(n=n, k=2)
f <- function(x) {sin(2*pi*x[1]) + .5*sin(4*pi*x[1]) +rnorm(1,0,.03) + x[2]^2}
y <- apply(x, 1, f) #f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian_beta$new(c(1, 1)), parallel=FALSE, verbose=10, nug.est=T)
ContourFunctions::cf(gp$predict, pts=x)
numDeriv::grad(func = function(x)gp$deviance(trend_params=x[1], params = x[2:4], nuglog=x[5]), x=c(-1.2, 1,1, 1, -4))
gp$deviance_grad(params = c(1,1,1), nug.update=T, nuglog=-4, trend_params=-1.2)
