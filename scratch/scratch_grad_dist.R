# 2D test
n <- 40
x <- matrix(runif(n*2), ncol=2)
f1 <- function(a) {sin(2*pi*a[1]) + sin(6*pi*a[2])}
#f1 <- TestFunctions::branin
#f1 <- TestFunctions::RFF_get(D=2)
y <- apply(x,1,f1) #+ rnorm(n,0,.01)
system.time(ContourFunctions::cf_data(x,y))
gp <- GauPro(x,y, verbose=2);gp$theta
system.time(ContourFunctions::cf_func(gp$pred, pts=x))
plot(y,gp$pred(x));abline(a=0,b=1)

xx <- runif(2)
gp$grad_dist(XX=xx)


# Check sample function
gp <- GauPro(x,y, verbose=2);gp$theta
xx <- runif(2)
samp1 <- gp$sample(xx, n=5000)
mean(samp1)
var(samp1)
gp$pred(xx, se.fit = T)
samp2 <- gp$sample(rbind(xx, c(.2,.3)), n=10000)
colMeans(samp2)
var(samp2)
gp$pred(rbind(xx, c(.2,.3)), se.fit = T)


# This checks that grad_dist matches up with numerical results
xx <- runif(2)
gp <- GauPro(x,y)
gp <- GauPro_kernel_model$new(x,y,kernel=Gaussian)
gp$grad_dist(xx)
gp$grad(xx)
eps <- 1e-4
xx3 <- rbind(xx, xx+c(eps,0), xx+c(0, eps))
# Get samples to estimate grad
samp3 <- gp$sample(XX = xx3, n = 1e5)
colMeans(samp3)
var(samp3)
# Estimate gradient in each dimension
grad_est1 <- (samp3[,2] - samp3[,1]) / eps
grad_est2 <- (samp3[,3] - samp3[,1]) / eps
summary(grad_est1)
gp$grad_dist(xx)
c(mean(grad_est1), mean(grad_est2)) #var(grad_est1)
var(cbind(grad_est1, grad_est2))

gp$grad_dist(matrix(c(.2,.3),ncol=2))
gp$grad_dist(matrix(c(.4,.5),ncol=2))
gp$grad_dist(matrix(c(.2,.3,.4,.5),byrow=T,ncol=2))


# > gp$grad_dist(matrix(c(.2,.3,.4,.5),byrow=T,ncol=2))
# $mean
# [,1]       [,2]
# [1,]  2.650846 -0.5852923
# [2,] -4.667038 -0.7475923
#
# $cov
# , , 1
#
# [,1]          [,2]
# [1,] 5.260113  0.0170742872
# [2,] 3.821696 -0.0002883557
#
# , , 2
#
# [,1]      [,2]
# [1,]  0.0170742872 0.2635080
# [2,] -0.0002883557 0.2533588

# Check grad dist and grad norm2 dist
tx <- matrix(runif(2),ncol=2)
# Next two should be equal to third
gp$grad_sample(tx, n=1e5) %>% colMeans()
gp$grad_sample(tx, n=1e5) %>% cov
gp$grad_dist(tx)
# Check grad_norm2_dist
gp$grad_sample(tx, n=1e5)^2 %>% rowSums %>% {c(mean(.), var(.))}
gp$grad_norm2_dist(tx)
sum(gp$grad(tx)^2) +gp$grad_dist(tx)$cov[1,,] %>% eigen() %>% .$val %>% sum
