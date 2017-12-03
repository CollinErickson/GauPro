# LOO test for kernel model
set.seed(0)
n <- 8
x <- matrix(seq(0,1,length.out = n), ncol=1)
f <- Vectorize(function(x) {sin(2*pi*x) + .5*sin(4*pi*x) +rnorm(1,0,.03)})
y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
#y[5] <- -.6
gp <- GauPro_Gauss_LOO$new(X=x, Z=y, verbose=10, nug.est=T)
gp <- GauPro_kernel_model_LOO$new(X=x, Z=y, verbose=10, nug.est=T, kernel=Matern52)
gp$plot1D()
gp$use_LOO <- F
gp$plot1D()



d <- 3
f <- TestFunctions::beambending
n <- 10*d
X <- (lhs::maximinLHS(n=n,k=d))
Z <- f(X)
nn <- 1e4
XX <- matrix(runif(nn*d),ncol=d)
ZZ <- f(XX)
gp1 <- GauPro_kernel_model$new(X=X, Z=Z, verbose=10, nug.est=T, kernel=Matern52)
gp2 <- GauPro_kernel_model_LOO$new(X=X, Z=Z, verbose=10, nug.est=T, kernel=Matern52)
p1 <- gp1$predict(XX, se=T)
p2 <- gp2$predict(XX, se=T)
plot(abs(ZZ-p1$mean), p1$se);abline(a=0,b=1,col=2)
plot(abs(ZZ-p2$mean), p2$se);abline(a=0,b=1,col=2)

tc <- comparer::mbc(gp1=GauPro_kernel_model$new(X=X, Z=Z, verbose=10, nug.est=T, kernel=Matern52),
              gp2=GauPro_kernel_model_LOO$new(X=X, Z=Z, verbose=10, nug.est=T, kernel=Matern52),
              inputi={X <- lhs::maximinLHS(n=n,k=d);Z <- f(X)},
              # targetin=list(XX,ZZ), target="ZZ"
              post=function(mod) {mod$predict(XX, se=T)}, target=ZZ,
              times=10, metric=c('rmse', 'mis90', 'sr27')
              )
print(tc)
plot(tc)
