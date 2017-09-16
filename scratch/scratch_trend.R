# trend_c
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

# trend_c 2D
set.seed(0)
n <- 60
x <- lhs::maximinLHS(n=n, k=2)
f <- TestFunctions::banana#Vectorize(function(x) {sin(2*pi*x) + .5*sin(4*pi*x) +rnorm(1,0,.3)}+10*x)
y <- 0*123 + f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian_beta$new(c(0,0)), parallel=FALSE, verbose=10, nug.est=T)
ContourFunctions::cf(gp$predict, batchmax=Inf, pts=x)
numDeriv::grad(func = function(x)gp$deviance(params=x[2:3], nuglog=x[4], trend_params=x[1]), x=c(10, 2,1, -4))
gp$deviance_grad(params = c(2,1), nug.update=T, nuglog=-4, trend_params=10)
numDeriv::grad(func = function(x)gp$deviance(trend_params=x[1], params=x[2:3], nuglog=x[4]), x=c(gp$trend$m, gp$kernel$beta, gp$kernel$logs2, log(gp$nug,10)))
gp$deviance_grad(params = c(gp$kernel$beta, gp$kernel$logs2), trend_params=gp$trend$m, nug.update=T, nuglog=log(gp$nug,10))
# Check fngr
gp$deviance(params = c(2,1), nuglog=-4, trend_params=10)
gp$deviance_grad(params = c(2,1), nug.update=T, nuglog=-4, trend_params=10)
gp$deviance_fngr(params = c(2,1), nug.update=T, nuglog=-4, trend_params=10)




# trend_LM
set.seed(0)
n <- 20
x <- matrix(seq(0,1,length.out = n), ncol=1)
f <- Vectorize(function(x) {sin(2*pi*x) + .5*sin(4*pi*x) +rnorm(1,0,.3)}+10*x)
y <- 123 + f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian$new(1), trend=trend_LM$new(D=1), parallel=FALSE, verbose=10, nug.est=T)
gp$cool1Dplot()
numDeriv::grad(func = function(x)gp$deviance(params=x[2:3], nuglog=x[4], trend_params=x[1]), x=c(10, 2,1, -4))
gp$deviance_grad(params = c(2,1), nug.update=T, nuglog=-4, trend_params=10)
numDeriv::grad(func = function(x)gp$deviance(trend_params=x[1], params=x[2:3], nuglog=x[4]), x=c(gp$trend$m, gp$kernel$beta, gp$kernel$logs2, log(gp$nug,10)))
gp$deviance_grad(params = c(gp$kernel$beta, gp$kernel$logs2), trend_params=gp$trend$m, nug.update=T, nuglog=log(gp$nug,10))
