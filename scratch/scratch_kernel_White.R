# Check numerically that gradient is correct for 1D
# For White kernel
set.seed(0)
n <- 20
x <- matrix(seq(0,1,length.out = n), ncol=1)
f <- Vectorize(function(x) {sin(2*pi*x) + .5*sin(4*pi*x) +rnorm(1,0,.3)})
y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=White$new(1), parallel=FALSE, verbose=10, nug.est=T)
gp$cool1Dplot()
numDeriv::grad(func = function(x)gp$deviance(params=x[2], trend_params=x[1], nuglog=x[3]), x=c(-.8, 2,-1))
gp$deviance_grad(params = c(2), nug.update=T, nuglog=-1, trend_params=-.8)
numDeriv::grad(func = function(x)gp$deviance(params=x[2], trend_params=x[1], nuglog=x[3]), x=c(gp$trend$m,  gp$kernel$logs2, log(gp$nug,10)))
gp$deviance_grad(params = c(gp$kernel$logs2), nug.update=T, nuglog=log(gp$nug,10), trend_params=gp$trend$b)


# Check if get same results as nugget

gp1 <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian$new(1), parallel=FALSE, verbose=10, nug.est=T)
gp2 <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian$new(1) + White$new(1), parallel=FALSE, verbose=10, nug.est=F)

c(gp1$kernel$beta, gp2$kernel$k1$beta)
c(gp1$nug * gp1$kernel$s2, gp2$nug * gp2$kernel$s2 + gp2$kernel$k2$s2)
# Looks like it works

gp3 <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian$new(1), parallel=FALSE, verbose=10, nug.est=F, nug=1e-8)
gp3$cool1Dplot()
gp2$cool1Dplot()





# Check dC_dtheta
m1 <- (gp$kernel$k(gp$X, beta=1,alpha=5) - gp$kernel$k(gp$X, beta=1-1e-6, alpha=5)) / 1e-6
C_nonug <- gp$kernel$k(gp$X, beta=1, alpha=5)
C <- C_nonug + gp$kernel$s2 * diag(gp$nug, nrow(C_nonug))
m2 <- gp$kernel$dC_dparams(params = c(1, 5, 1), X = gp$X, C = C, C_nonug = C_nonug)[1,,]
summary(c(m1-m2))# %>% summary
plot(m1, m2)

gp$deviance_grad()
dsign <- 1
mm1 <- gp$kernel$dC_dparams(C = C, C_nonug = C_nonug, X=gp$X)[[1]]
dsign <- -1
mm2 <- gp$kernel$dC_dparams(C = C, C_nonug = C_nonug, X=gp$X)[[1]]
plot(c(mm1[[1]]), c(mm2[[1]]))


# Check dC_dlogs2
beta <- gp$kernel$beta
s2 <- gp$kernel$s2+.2
nug <- gp$nug
eps <- 1e-6
m1 <- (gp$kernel$k(gp$X, beta=beta, s2=s2+eps) - gp$kernel$k(gp$X, beta=beta, s2=s2-eps)) / eps / 2
C_nonug <- gp$kernel$k(gp$X, beta=beta, s2=s2)
C <- C_nonug + s2 * diag(nug, nrow(C_nonug))
m2 <- gp$kernel$dC_dparams(params = log(c(beta, s2),10), X = gp$X, C = C, C_nonug = C)[[1]][[2]]
c(m1 * s2 * log(10) -m2) %>% summary
plot(c(m1 * s2 * log(10)), c(m2))


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
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=RatQuad$new(c(1, 1)), parallel=FALSE, verbose=10, nug.est=T)
ContourFunctions::cf(gp$predict, pts=x, batchmax=Inf)
ContourFunctions::cf(f, pts=x)
numDeriv::grad(func = function(x)gp$deviance(params = x[1:3], nuglog=x[4]), x=c(1,1, 1, -4))
gp$deviance_grad(params = c(1,1,1), nug.update=T, nuglog=-4)
