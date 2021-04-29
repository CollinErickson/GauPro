n <- 5
for (i in 0:(n-1-1)) {
  for (j in (i+1):(n-1)) {
    k = (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1
    print(c(i, j, k))
  }
}

kk <- IndexKernel$new(D=1, nlevels=5, xindex=1)
kk
kk$kone(1,1, s2=1)
# kk$kone(1,2, s2=1)
kk$k(1,1)
kk$k(1,2)

kk$p <- (1:10)/100
kmat <- outer(1:5, 1:5, Vectorize(kk$k))
kmat <- matrix(NA,5,5)
for (i in 1:5) {
  for (j in 1:5) {
    kmat[i,j] <- kk$k(i,j)
    #print(c(i,j, kk$k(i,j)))
  }
}
kmat

kk$k(matrix(c(1:5),ncol=1))
kk$p[2] <- .1
kk$k(matrix(c(1:5),ncol=1))

kk$param_optim_lower()
kk$param_optim_upper()
kk$param_optim_start()
kk$param_optim_start0()



# 1D example
library(dplyr)
X <- matrix(sample(1:3, size=12, replace=T), ncol=1)
Z <- c(X) %>% {ifelse(.==2,10,0) + rnorm(length(.),0,.2)}
tibble(X=X[,1], Z) %>% arrange(X,Z)
k2 <- IndexKernel$new(D=1, nlevels=3, xind=1)
# debugonce(k2$dC_dparams)
# k2$p_upper <- .99*k2$p_upper
gp <- GauPro_kernel_model$new(X=X, Z=Z, kernel = k2, verbose = 5)
gp$kernel$p
gp$kernel$k(x = gp$X)
tibble(X=X[,1], Z=Z, pred=gp$predict(X)[,1]) %>% arrange(X, Z)
tibble(X=X[,1], Z) %>% group_by(X) %>% summarize(n=n(), mean(Z))
# gp$plot1D(xmin = 1.5, xmax=2.5)


# 2D, kernel on 1D
n <- 20
X <- cbind(matrix(sample(1:3, size=n, replace=T), ncol=1),
           matrix(sample(1:3, size=n, replace=T), ncol=1))
Z <- X[,1] - X[,2] + rnorm(n,0,.1)
tibble(X=X, Z) %>% arrange(X,Z)
k2 <- IndexKernel$new(D=2, nlevels=3, xind=1)
# debugonce(k2$dC_dparams)
# k2$p_upper <- .65*k2$p_upper
gp <- GauPro_kernel_model$new(X=X, Z=Z, kernel = k2, verbose = 5)
gp$kernel$p
gp$kernel$k(x = gp$X)
tibble(X=X, Z=Z, pred=gp$predict(X)[,1]) %>% arrange(X, Z)
tibble(X=X[,1], Z) %>% group_by(X) %>% summarize(n=n(), mean(Z))


# 2D, kernel on both dim
n <- 12
X <- cbind(matrix(sample(1:3, size=n, replace=T), ncol=1),
           matrix(sample(1:3, size=n, replace=T), ncol=1))
Z <- X[,1] - X[,2]^2 + rnorm(n,0,.1)
tibble(X=X, Z) %>% arrange(X,Z)
k2 <- IndexKernel$new(D=2, nlevels=3, xind=1) * IndexKernel$new(D=2, nlevels=3, xind=2)
# debugonce(k2$dC_dparams)
k2$k1$p_upper <- .5*k2$k1$p_upper
k2$k2$p_upper <- .5*k2$k2$p_upper
gp <- GauPro_kernel_model$new(X=X, Z=Z, kernel = k2, verbose = 5, nug.min = 1e-4)
gp$kernel$k1$p
gp$kernel$k2$p
gp$kernel$k(x = gp$X)
tibble(X=X, Z=Z, pred=gp$predict(X)[,1]) %>% arrange(X, Z)
tibble(X=X, Z) %>% group_by(X) %>% summarize(n=n(), mean(Z))
XX <- as.matrix(expand.grid(1:3,1:3))
bind_cols(XX, gp$pred(XX))
cbind(XX, gp$kernel$k(x = XX))



# 2D, Gaussian on 1D, index on 2nd dim
library(dplyr)
n <- 20
X <- cbind(matrix(runif(n,2,6), ncol=1),
           matrix(sample(1:2, size=n, replace=T), ncol=1))
X <- rbind(X, c(3.3,3))
n <- nrow(X)
Z <- X[,1] - (X[,2]-1.8)^2 + rnorm(n,0,.1)
tibble(X=X, Z) %>% arrange(X,Z)
k2a <- IgnoreIndsKernel$new(k=Gaussian$new(D=1), ignoreinds = 2)
k2b <- IndexKernel$new(D=2, nlevels=3, xind=2)
k2 <- k2a * k2b
# debugonce(k2$dC_dparams)
k2b$p_upper <- .65*k2b$p_upper
gp <- GauPro_kernel_model$new(X=X, Z=Z, kernel = k2, verbose = 5, nug.min=1e-2)
gp$kernel$k1$kernel$beta
gp$kernel$k2$p
gp$kernel$k(x = gp$X)
tibble(X=X, Z=Z, pred=gp$predict(X)[,1]) %>% arrange(X, Z)
tibble(X=X[,2], Z) %>% group_by(X) %>% summarize(n=n(), mean(Z))
curve(gp$pred(cbind(matrix(x,ncol=1),1)),2,6, ylim=c(min(Z), max(Z))); points(X[X[,2]==1,1], Z[X[,2]==1])
curve(gp$pred(cbind(matrix(x,ncol=1),2)), add=T, col=2); points(X[X[,2]==2,1], Z[X[,2]==2], col=2)
curve(gp$pred(cbind(matrix(x,ncol=1),3)), add=T, col=3); points(X[X[,2]==3,1], Z[X[,2]==3], col=3)
legend(legend=1:3, fill=1:3, x="topleft")
cbind(X, cov=gp$kernel$k(X, c(5.5,3))) %>% arrange(-cov)
# See which points affect (5.5, 3 themost)
data.frame(X, cov=gp$kernel$k(X, c(5.5,3))) %>% arrange(-cov)
