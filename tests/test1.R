source("R/corr.R")
gp <- GauPro$new()
gp$fit(matrix(runif(6),3,2),1:3)
gp$pred(matrix(runif(6),3,2))

# 1D test
n <- 12
x <- matrix(seq(0,1,length.out = 12), ncol=1)
y <- sin(2*pi*x) #+ rnorm(n,0,1e-1)
y <- sqrt(x)-x
y <- (2*x) %%1
y <- c(y)
plot(x,y)
gp <- GauPro$new()
gp$fit(x,y)
#gp$cool1Dplot()
#gp$pred(matrix(runif(6),6,1))
curve(gp$pred(x)$me);points(x,y)
curve(gp$pred(x)$me+2*gp$pred(x)$se,col=2,add=T);curve(gp$pred(x)$me-2*gp$pred(x)$se,col=2,add=T)
#gpf <- UGP::UGP$new(X=x,Z=y, package="mlegp")
#gpf$mod[[1]]$mu;gpf$mod[[1]]$beta;gpf$mod[[1]]$nug
#gp$mu_hat
# deviance test
curve(sapply(x, gp$deviance_log),-10,10, n = 300)
gp$deviance_search()

gp$theta
gp$deviance_search()
gp$theta_update()
gp$theta

gp$deviance_search3()
c(gp$theta,gp$nug)
gp$all_update()
c(gp$theta,gp$nug)

gp <- GauPro$new();gp$fit(x,y);gp$cool1Dplot()

# gradtest
gp <- GauPro$new()
gp$fit(x,y)
gp$mean_grad(c(.5,.6))
par(mfrow=c(2,1))
curve(gp$pred(x)$me);points(x,y)
curve(gp$pred(x)$me+2*gp$pred(x)$se,col=2,add=T);curve(gp$pred(x)$me-2*gp$pred(x)$se,col=2,add=T)
curve(gp$mean_grad(x), add=F, col=3)
abline(v=x)
curve(gp$mean_grad_norm(x))
gp$mean_grad_norm(c(.5,.6))

# 2D test
n <- 40
x <- matrix(runif(n*2), ncol=2)
f1 <- function(a) {sin(3*pi*a[1]) + sin(3*pi*a[2])}
y <- apply(x,1,f1)
system.time(contourfilled::contourfilled.data(x,y))
gp <- GauPro$new()
gp$fit(x,y);gp$theta
gp$pred(matrix(runif(6),3,2))
c(gp$mu_hat, gp$s2_hat)
system.time(contourfilled::contourfilled.func(gp$pred,out.name = 'mean', pts=x))
plot(y,gp$pred(x)$mean);abline(a=0,b=1)
system.time(print(gp$deviance_search()))
system.time(print(gp$deviance_search3()))
c(gp$theta, gp$nug)
gp$all_update()
c(gp$theta, gp$nug)
system.time(contourfilled::contourfilled.func(gp$pred,out.name = 'mean', pts=x))

gp$mean_grad(matrix(c(.5, .75, .5, .83),2,2,byrow=T))
gp$mean_grad(matrix(c(.5, .75),1,2,byrow=T))
contourfilled::contourfilled.func(function(xx){gp$mean_grad(xx)[2]})
gp$mean_grad_norm(matrix(c(.5, .75, .5, .83),2,2,byrow=T))
gp$mean_grad_norm(c(.5,.5))
contourfilled::contourfilled.func(function(xx){gp$mean_grad_norm(xx)})


gpf <- UGP::UGP$new(X=x,Z=y, package="mlegp")
gpf$mod[[1]]




# higher dim test
n <- 60
d <- 3
x <- matrix(runif(n*d), ncol=d)
f1 <- function(a) {sum(sin(1:d*pi*a))}
y <- apply(x,1,f1)
gp <- GauPro$new()
gp$fit(x,y);c(gp$theta,gp$nug)
nn <- 20
gp$pred(matrix(runif(nn*d),ncol=d))
gp$mean_grad(matrix(runif(nn*d),ncol=d))
gp$mean_grad_norm(matrix(runif(nn*d),ncol=d))
plot(y,gp$pred(x)$mean);abline(a=0,b=1)
