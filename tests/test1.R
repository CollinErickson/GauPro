#source("R/corr.R")
#gp <- GauPro$new()
#gp$fit(matrix(runif(6),3,2),1:3)
#gp$pred(matrix(runif(6),3,2))

# 1D test
n <- 12
x <- matrix(seq(0,1,length.out = 12), ncol=1)
y <- sin(2*pi*x) #+ rnorm(n,0,1e-1)
y <- sqrt(x)-x
y <- (2*x) %%1
y <- c(y)
plot(x,y)
gp <- GauPro$new(X=x, Z=y)
curve(gp$pred(x));points(x,y)
curve(gp$pred(x)+2*gp$pred(x,T)$se,col=2,add=T);curve(gp$pred(x)-2*gp$pred(x,T)$se,col=2,add=T)
curve(sapply(x, gp$deviance_theta_log),-10,10, n = 300) # deviance profile
gp$optim()

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




# 2D test
n <- 40
x <- matrix(runif(n*2), ncol=2)
f1 <- function(a) {sin(3*pi*a[1]) + sin(3*pi*a[2])}
y <- apply(x,1,f1)
system.time(contourfilled::contourfilled.data(x,y))
gp <- GauPro$new(x,y);gp$theta
gp$pred(matrix(runif(6),3,2))
c(gp$mu_hat, gp$s2_hat)
system.time(contourfilled::contourfilled.func(gp$pred, pts=x))
plot(y,gp$pred(x));abline(a=0,b=1)
#system.time(print(gp$deviance_search()))
#system.time(print(gp$deviance_search3()))
c(gp$theta, gp$nug)
gp$update()
c(gp$theta, gp$nug)
system.time(contourfilled::contourfilled.func(gp$pred, pts=x))

gp$grad(matrix(c(.5, .75, .5, .83),2,2,byrow=T))
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
n <- 60
d <- 3
x <- matrix(runif(n*d), ncol=d)
f1 <- function(a) {sum(sin(1:d*pi/a) + (1/a))}
y <- apply(x,1,f1)
gp <- GauPro$new(x,y, verbose=2);c(gp$theta,gp$nug)
nn <- 20
gp$pred(matrix(runif(nn*d),ncol=d))
gp$grad(matrix(runif(nn*d),ncol=d))
gp$grad_norm(matrix(runif(nn*d),ncol=d))
plot(y,gp$pred(x));abline(a=0,b=1)
