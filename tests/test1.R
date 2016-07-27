source("R/corr.R")
gp <- GauPro$new()
gp$fit(matrix(runif(6),3,2),1:3)
gp$pred(matrix(runif(6),3,2))

# 1D test
n <- 12
x <- matrix(seq(0,1,length.out = 12), ncol=1)
y <- sin(2*pi*x) + rnorm(n,0,1e-1)
y <- sqrt(x)-x
y <- (2*x) %%1
y <- c(y)
plot(x,y)
gp <- GauPro$new()
gp$fit(x,y)
gp$cool1Dplot()
gp$pred(matrix(runif(6),6,1))
curve(gp$pred(x)$me);points(x,y)
curve(gp$pred(x)$me+2*gp$pred(x)$se,col=2,add=T);curve(gp$pred(x)$me-2*gp$pred(x)$se,col=2,add=T)
#gpf <- UGP::UGP$new(X=x,Z=y, package="mlegp")
#gpf$mod[[1]]$mu;gpf$mod[[1]]$beta;gpf$mod[[1]]$nug
#gp$mu_hat
# deviance test
curve(sapply(x, gp$deviance_log),-10,10)
gp$deviance_search()

gp$theta
gp$deviance_search()
gp$theta_update()
gp$theta

gp <- GauPro$new();gp$fit(x,y);gp$cool1Dplot()

# 2D test
n <- 40
x <- matrix(runif(n*2), ncol=2)
f1 <- function(a) {sin(1*pi*a[1]) + sin(6*pi*a[2])}
y <- apply(x,1,f1)
system.time(contourfilled::contourfilled.data(x,y))
gp <- GauPro$new()
gp$fit(x,y)
gp$pred(matrix(runif(6),3,2))
rbind(gp$mu_hat, gp$s2_hat)
system.time(contourfilled::contourfilled.func(gp$pred,out.name = 'mean', pts=x))
plot(y,gp$pred(x)$mean);abline(a=0,b=1)
gpf <- UGP::UGP$new(X=x,Z=y, package="mlegp")
gpf$mod[[1]]
