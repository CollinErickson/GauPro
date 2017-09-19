# Making sure that making covariance functions into Rcpp functions
#  still gives correct results.

# Exponential kernel

set.seed(0)
n <- 20
x <- matrix(seq(0,1,length.out = n), ncol=1)
f <- Vectorize(function(x) {sin(2*pi*x) + .5*sin(4*pi*x) +rnorm(1,0,.3)})
y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
system.time(gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Exponential$new(1), parallel=FALSE, verbose=10, nug.est=T))
# .89 sec
system.time(gp$cool1Dplot()) # .42 sec


gp$predict(.656) # -0.6040612
gp$predict(c(.11, .24, .455, .676, .888)) # 1.5120375, 0.8360396, 0.4850529, -0.6252635, -1.3454632
gp$predict(matrix(c(.11, .24, .455, .676, .888), ncol=1))

set.seed(0)
n <- 200
x <- matrix(runif(6*n), ncol=6)
y <- TestFunctions::OTL_Circuit(x)
system.time(gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Exponential, parallel=FALSE, verbose=10, nug.est=T))
#  19.68 / 20.28 s
system.time(gp$predict(x+.01)) # .43 sec
system.time(gp$predict(x+.01, covmat = T)) # .72 sec
gp$predict(matrix(c(.1,.2,.3,.4,.5,.6), ncol=6)) # 5.577286







# Matern 3/2 kernel

set.seed(0)
n <- 20
x <- matrix(seq(0,1,length.out = n), ncol=1)
f <- Vectorize(function(x) {sin(2*pi*x) + .5*sin(4*pi*x) +rnorm(1,0,.3)})
y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
system.time(gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Matern32, parallel=FALSE, verbose=10, nug.est=T))
# 1.73 sec
system.time(gp$cool1Dplot()) # .55 sec


gp$predict(.656) # -0.6063402
gp$predict(c(.11, .24, .455, .676, .888)) # 1.4436862  0.8492838  0.4596046 -0.6550763 -1.2473287
gp$predict(matrix(c(.11, .24, .455, .676, .888), ncol=1))

set.seed(0)
n <- 200
x <- matrix(runif(6*n), ncol=6)
y <- TestFunctions::OTL_Circuit(x)
system.time(gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Matern32, parallel=FALSE, verbose=10, nug.est=T))
#  29.31 / 30.49 s
system.time(gp$predict(x+.01)) # .65 sec
system.time(gp$predict(x+.01, covmat = T)) # 1.15 sec
gp$predict(matrix(c(.1,.2,.3,.4,.5,.6), ncol=6)) # 5.646576


# Matern 5/2 kernel

set.seed(0)
n <- 20
x <- matrix(seq(0,1,length.out = n), ncol=1)
f <- Vectorize(function(x) {sin(2*pi*x) + .5*sin(4*pi*x) +rnorm(1,0,.3)})
y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
system.time(gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Matern52, parallel=FALSE, verbose=10, nug.est=T))
# 1.59 sec
system.time(gp$cool1Dplot()) # .56 sec


gp$predict(.656) # -0.616631
gp$predict(c(.11, .24, .455, .676, .888)) # 1.4023642 0.8733849 0.4285692 -0.6816842-1.1858629
gp$predict(matrix(c(.11, .24, .455, .676, .888), ncol=1))

set.seed(0)
n <- 200
x <- matrix(runif(6*n), ncol=6)
y <- TestFunctions::OTL_Circuit(x)
system.time(gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Matern52, parallel=FALSE, verbose=10, nug.est=T))
#  24.51 / 25.66 s
system.time(gp$predict(x+.01)) # .68 sec
system.time(gp$predict(x+.01, covmat = T)) # 1.02 sec
gp$predict(matrix(c(.1,.2,.3,.4,.5,.6), ncol=6)) # 5.526564
