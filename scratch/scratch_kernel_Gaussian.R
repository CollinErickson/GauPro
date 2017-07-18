# Check numerically that gradient is correct for 1D
set.seed(0)
n <- 20
x <- matrix(seq(0,1,length.out = n), ncol=1)
f <- Vectorize(function(x) {sin(2*pi*x) + .5*sin(4*pi*x) +rnorm(1,0,.03)})
y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian$new(1), parallel=FALSE, verbose=10, nug.est=T)
gp$cool1Dplot()
numDeriv::grad(func = function(x)gp$deviance(params = x[1:2], nuglog=x[3]), x=c(100,1, -4))
gp$deviance_grad(params = c(100,1), nug.update=T, nuglog=-4)
numDeriv::grad(func = function(x)gp$deviance(params=x[1:2], nuglog=x[3]), x=c(gp$kernel$theta, gp$kernel$s2, log(gp$nug,10)))
gp$deviance_grad(params = c(gp$kernel$theta, gp$kernel$s2), nug.update=T, nuglog=log(gp$nug,10))

# Check dC_dtheta
m1 <- (gp$kernel$k(gp$X, theta=100) - gp$kernel$k(gp$X, theta=100-1e-6)) / 1e-6
C <- gp$kernel$k(gp$X, theta=100)
m2 <- gp$kernel$dC_dparams(params = c(100, 1), X = gp$X, C = C, C_nonug = C)[[1]][[1]]
c(m1-m2) %>% summary

# Check if optim works
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian$new(1), parallel=FALSE, verbose=10)
# gp$optim()
rmse <- function(gp) {
  xx <- seq(0,1,l=201)
  zz <- f(xx)
  sqrt(mean((gp$predict(xx) - zz)^2))
}
rmse(gp)

# Check nug with grad
numDeriv::grad(func = function(x) {gp$deviance(nuglog = x[3], params = x[1:2])}, x=c(100,1, -6))
gp$deviance_grad(params = c(100,1), nug.update=T, nuglog = -6)





# Check 2D
set.seed(0)
n <- 30
x <- lhs::maximinLHS(n=n, k=2)
f <- function(x) {sin(2*pi*x[1]) + .5*sin(4*pi*x[1]) +rnorm(1,0,.03) + x[2]^2}
y <- apply(x, 1, f) #f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian$new(c(1, 1)), parallel=FALSE, verbose=10, nug.est=T)
ContourFunctions::cf(gp$predict)
numDeriv::grad(func = gp$deviance, x=c(100,1))
gp$deviance_grad(params = c(100,1), nug.update=F)
