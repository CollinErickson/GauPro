# Check numerically that gradient is correct for 1D
# For Gaussian_logl kernel
set.seed(0)
n <- 20
x <- matrix(seq(0,1,length.out = n), ncol=1)
f <- Vectorize(function(x) {sin(2*pi*x) + .5*sin(4*pi*x) +rnorm(1,0,.3)})
y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian_logl$new(1), parallel=FALSE, verbose=10, nug.est=T)
gp$cool1Dplot()
numDeriv::grad(func = function(x)gp$deviance(params=x[1:2], nuglog=x[3]), x=c(2,1, -4))
gp$deviance_grad(params = c(2,1), nug.update=T, nuglog=-4)
numDeriv::grad(func = function(x)gp$deviance(params=x[1:2], nuglog=x[3]), x=c(gp$kernel$logl, gp$kernel$logs2, log(gp$nug,10)))
gp$deviance_grad(params = c(gp$kernel$logl, gp$kernel$logs2), nug.update=T, nuglog=log(gp$nug,10))

# Check dC_dtheta
m1 <- (gp$kernel$k(gp$X, logl=1) - gp$kernel$k(gp$X, logl=1-1e-6)) / 1e-6
C <- gp$kernel$k(gp$X, logl=1)
m2 <- gp$kernel$dC_dparams(params = c(1, 1), X = gp$X, C = C, C_nonug = C)[[1]][[1]]
c(m1-m2) %>% summary




# Check 2D
set.seed(0)
n <- 30
x <- lhs::maximinLHS(n=n, k=2)
f <- function(x) {sin(2*pi*x[1]) + .5*sin(4*pi*x[1]) +rnorm(1,0,.03) + x[2]^2}
y <- apply(x, 1, f) #f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian_beta$new(c(1, 1)), parallel=FALSE, verbose=10, nug.est=T)
ContourFunctions::cf(gp$predict, pts=x)
numDeriv::grad(func = function(x)gp$deviance(params = x[1:3], nuglog=x[4]), x=c(1,1, 1, -4))
gp$deviance_grad(params = c(1,1,1), nug.update=T, nuglog=-4)
