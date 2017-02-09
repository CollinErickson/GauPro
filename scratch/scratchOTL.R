

# higher dim test
n <- 21
d <- 6
x <- matrix(runif(n*d), ncol=d)
f1 <- function(a) {sum(sin(1:d*pi/a))}
f1 <- function(a) {sum(sin(1*pi/a[1:5]))}
f1 <- TestFunctions::OTL_Circuit
y <- apply(x,1,f1) #+ rnorm(n,0,.1)
system.time(gp <- GauPro(x,y, verbose=0));c(gp$theta,gp$nug)
summary(gp$predict(matrix(runif(6*1e3), ncol=6)))

microbenchmark::microbenchmark(GauPro$new(x,y, useOptim2=F), GauPro$new(x,y, useOptim2=T), times = 10)
microbenchmark::microbenchmark(GauPro$new(x,y), times = 100)
nn <- 2000
gp$pred(matrix(runif(nn*d),ncol=d))
gp$grad(matrix(runif(nn*d),ncol=d))
gp$grad_norm(matrix(runif(nn*d),ncol=d))
plot(y,gp$pred(x));abline(a=0,b=1)

mod <- UGP::IGP(X=x,Z=y, package='GauPro')
summary(mod$predict(matrix(runif(6*1e3), ncol=6)))
