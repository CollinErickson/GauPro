# Making sure that making covariance functions into Rcpp functions
#  still gives correct results.

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
