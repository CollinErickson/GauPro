# 1D test
n <- 12
x <- matrix(seq(0,1,length.out = n), ncol=1)
y <- sin(2*pi*x) + rnorm(n,0,1e-1)
y <- sqrt(x)-x
y <- (2*x) %%1
plot(x,y)
gp <- GauPro$new(X=x, Z=y, useLLH=F, useGrad=F)
curve(gp$pred(x));points(x,y)
curve(gp$pred(x)+2*gp$pred(x,T)$se,col=2,add=T);curve(gp$pred(x)-2*gp$pred(x,T)$se,col=2,add=T)
curve(sapply(x, gp$deviance_theta_log),-10,10, n = 300) # deviance profile
gp$optim()
microbenchmark::microbenchmark(GauPro$new(x,y, useOptim2=F), GauPro$new(x,y, useOptim2=T), times = 100)

gp$optim()
c(gp$theta,gp$nug)
gp$update()
c(gp$theta,gp$nug)

gp <- GauPro$new(x,y);gp$cool1Dplot()

# gradtest
gp <- GauPro$new(x,y)
gp$grad(c(.5,.6))
par(mfrow=c(2,1))
curve(gp$pred(x));points(x,y)
curve(gp$pred(x)+2*gp$pred(x,T)$se,col=2,add=T);curve(gp$pred(x)-2*gp$pred(x,T)$se,col=2,add=T)
curve(gp$grad(x), add=F, col=3)
abline(v=x)
curve(gp$grad_norm(x))
gp$grad_norm(c(.5,.6))

#compare grad times
ngrad <- 1
dgrad <- ncol(x)
xgrad <- matrix(runif(ngrad*dgrad),ncol=dgrad)
system.time(gp$grad(xgrad))
# num no longer available system.time(gp$grad_num(xgrad))
#                     plot(gp$grad_norm(xgrad),gp$grad_num_norm(xgrad))

# Test update function
xnew <- matrix(runif(4),ncol=1)
znew <- (2*xnew) %%1
gp$cool1Dplot()
gp$update(Xnew=xnew,Znew=znew)
gp$cool1Dplot()

# time
library(lineprof)
l <- lineprof(GauPro$new(X=x,Z=y))
shine(l)
microbenchmark::microbenchmark(gp$optim(restarts = 32), gp$optimParallel(restarts = 32), times = 10)
microbenchmark::microbenchmark(gp$optim(), gp$optimParallel(), times = 10)
# compare useC and parallel options
gp1 <- GauPro$new(x,y,parallel=F,useC=T)
gp2 <- GauPro$new(x,y,parallel=F,useC=F)
gp3 <- GauPro$new(x,y,parallel=T,useC=F)
gp4 <- GauPro$new(x,y,parallel=T,useC=T)
microbenchmark::microbenchmark(gp1$optim(32), gp2$optim(32), gp3$optim(32), gp4$optim(32), times = 5)
gp1 <- GauPro$new(x,y,parallel=F,useC=T)
gp2 <- GauPro$new(x,y,parallel=T,useC=F)
gp3 <- GauPro$new(x,y,parallel=F,useC=F)

# Check deviance grad
curve(sapply(x, gp$deviance),1,10, n = 300,lwd=5) # deviance profile
curve(sapply(x, gp$deviance_grad),1,10, n = 300,lwd=5) # deviance profile
curve(sapply(x, function(xx)numDeriv::grad(gp$deviance,xx)),1,10, n = 300,col=2,add=T,lwd=3) # deviance profile
tx <- seq(1,10,length.out = 200)
plot(sapply(tx, gp$deviance_grad), sapply(tx, function(xx)numDeriv::grad(gp$deviance,xx)))
# on beta scale
curve(sapply(x, gp$deviance_log_grad),1,10, n = 300,lwd=5) # deviance profile
curve(sapply(x, function(xx)numDeriv::grad(gp$deviance_log,xx)),1,10, n = 300,col=2,add=T,lwd=3) # deviance profile
# time compare
microbenchmark::microbenchmark(gp$deviance_log_grad(2), numDeriv::grad(gp$deviance_log,2),times=1e4)
microbenchmark::microbenchmark(gp$deviance_log_grad(c(1,2)), numDeriv::grad(gp$deviance_log,c(1,2)),times=1e2)
microbenchmark::microbenchmark(GauPro$new(x,y, useOptim2=T,parallel=F, useGrad=F, verbose=0),
                               GauPro$new(x,y, useOptim2=T,parallel=F, useGrad=T, verbose=0),times=100)
# Check for nugget
nugdev <- Vectorize(function(xxx)gp$deviance(nug=xxx))
nuggrad <- Vectorize(function(xxx)gp$deviance_gradC(nug=xxx,overwhat="nug"))
nugnumgrad <- Vectorize(function(xx)numDeriv::grad(nugdev, xx))
curve(nugdev,1e-6,1e-1,log='x',n = 200)
curve(nuggrad,1e-6,1e-1,log='x',n = 200)
curve(nugnumgrad,1e-6,1e-1,log=c('x'),n=200)
xnug <- 10 ^ (seq(-8,-1,length.out = 1000))
plot(xnug[1:999], diff(nugdev(xnug))/diff(xnug), log=c("x","y"))
plot(diff(nugdev(xnug))/diff(xnug), nuggrad(xnug[1:999]))

# Check fngr
microbenchmark::microbenchmark(gp$deviance_fngr(), {gp$deviance();gp$deviance_gradC()}, {gp$devianceC();gp$deviance_gradC()})

# Check LLH and grad
curve(sapply(10^x, gp$deviance_LLH_grad),-10,10, n = 300,ylim=c(-4e7,10)) # deviance profile
curve(sapply(10^x, function(a)numDeriv::grad(gp$deviance_LLH,a)),-10,10, n = 300) # deviance profile
tx <- seq(-4,4,length.out = 200)
plot(sapply(10^tx, gp$deviance_LLH_grad), sapply(10^tx, function(a)numDeriv::grad(gp$deviance_LLH,a)))
numDeriv::grad(function(a)gp$deviance_LLH_log2(joint=a), log(c(gp$theta, gp$nug),10))
gp$deviance_LLH_log2_grad()

# Check deviance speeds
microbenchmark::microbenchmark(gp$deviance(), gp$devianceC(), gp$devianceCC(),times=1e4)


# 2D test
n <- 80
x <- matrix(runif(n*2), ncol=2)
f1 <- function(a) {sin(3*pi*a[1]) + sin(3*pi*a[2])}
#f1 <- TestFunctions::branin
#f1 <- TestFunctions::RFF_get(D=2)
y <- apply(x,1,f1) + rnorm(n,0,.01)
system.time(contourfilled::contourfilled.data(x,y))
gp <- GauPro$new(x,y, verbose=2, useLLH=T);gp$theta
system.time(contourfilled::contourfilled.func(gp$pred, pts=x))
plot(y,gp$pred(x));abline(a=0,b=1)
microbenchmark(GauPro$new(x,y), GauPro$new(x,y), times = 10)
#system.time(print(gp$deviance_search()))
#system.time(print(gp$deviance_search3()))
c(gp$theta, gp$nug)
gp$update()
c(gp$theta, gp$nug)
system.time(contourfilled::contourfilled.func(gp$pred, pts=x))

p$grad(matrix(c(.5, .75, .5, .83),2,2,byrow=T))
gp$grad(matrix(c(.5, .75),1,2,byrow=T))
contourfilled::contourfilled.func(function(xx){gp$grad(xx)[2]})
gp$grad_norm(matrix(c(.5, .75, .5, .83),2,2,byrow=T))
gp$grad_norm(c(.5,.5))
contourfilled::contourfilled.func(function(xx){gp$grad_norm(xx)})

#optim check, show contour of theta
contourfilled::contourfilled.func(function(bet) {gp$deviance_log(beta = bet)}, xcontlim = c(-10,10), ycontlim = c(-10,10))

gpf <- UGP::UGP$new(X=x,Z=y, package="mlegp")
gpf$mod[[1]]




# higher dim test
n <- 200
d <- 4
x <- matrix(runif(n*d), ncol=d)
f1 <- function(a) {sum(sin(1:d*pi/a))}
y <- apply(x,1,f1) + rnorm(n,0,.1)
gp <- GauPro$new(x,y, verbose=0, parallel=T, useC=F);c(gp$theta,gp$nug)
microbenchmark::microbenchmark(GauPro$new(x,y, useOptim2=F), GauPro$new(x,y, useOptim2=T), times = 10)
microbenchmark::microbenchmark(GauPro$new(x,y), times = 100)
nn <- 2000
gp$pred(matrix(runif(nn*d),ncol=d))
gp$grad(matrix(runif(nn*d),ncol=d))
gp$grad_norm(matrix(runif(nn*d),ncol=d))
plot(y,gp$pred(x));abline(a=0,b=1)
