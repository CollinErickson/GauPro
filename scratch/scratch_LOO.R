# LOO test for kernel model
set.seed(0)
n <- 8
x <- matrix(seq(0,1,length.out = n), ncol=1)
f <- Vectorize(function(x) {sin(2*pi*x) + .5*sin(4*pi*x) +rnorm(1,0,.03)})
y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
#y[5] <- -.6
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=RatQuad$new(1, 1.5), parallel=FALSE, verbose=10, nug.est=T)
gp$cool1Dplot()

# Plot LOOs on top of coolplot
loo <- gp$pred_LOO()
points(x, loo, col=2)
plot(y, loo);abline(a=0,b=1)

# Check to see if LOOs match predictions when row is actually removed
gp2store <- matrix(NA, n, 4)
for (i in 1:n) {
  gp2 <- gp$clone(deep = TRUE)
  gp2$update(Xall = x[-i, , drop=FALSE], Zall = y[-i], no_update = T)
  # gp2$cool1Dplot()
  gp2store[i, ] <- c(y[i], gp2$pred(x[i]), gp$pred(x[i]), loo[i]) #%>% print
}
pairs(gp2store)
gp2store
plot(gp2store[,2], gp2store[,4]);abline(a=0,b=1)
gp$cool1Dplot()
points(x, gp2store[,4], col=2)
points(x, gp2store[,2], col=3)

# Plot with se
gp$cool1Dplot()
loo <- gp$pred_LOO(se.fit=T)
points(x, loo$fit, col=2)
points(x, loo$fit + 2*loo$se, col=3)
points(x, loo$fit - 2*loo$se, col=3)






# LOO test for non-kernel model
set.seed(0)
n <- 8
x <- matrix(seq(0,1,length.out = n), ncol=1)
f <- Vectorize(function(x) {sin(2*pi*x) + .5*sin(4*pi*x) +rnorm(1,0,.03)})
y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
#y[5] <- -.6
gp <- GauPro_Gauss$new(X=x, Z=y, parallel=FALSE, verbose=10, nug.est=T)
gp$cool1Dplot()

# Plot LOOs on top of coolplot
loo <- gp$pred_LOO()
points(x, loo, col=2)
plot(y, loo);abline(a=0,b=1)

# Plot with se
gp$cool1Dplot()
loo <- gp$pred_LOO(se.fit=T)
points(x, loo$fit, col=2)
points(x, loo$fit + 2*loo$se, col=3)
points(x, loo$fit - 2*loo$se, col=3)
