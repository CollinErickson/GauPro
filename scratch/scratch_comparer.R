library(comparer)
n <- 30
d <- 2
f <- TestFunctions::banana
x <- matrix(runif(n*d), n,d)
y <- f(x)
nn <- 300
xx <- matrix(runif(nn*d), nn,d)
yy <- f(xx)
mbc(Gaussian, Matern32, Matern52, evaluator={GauPro_kernel_model$new(X=x,Z=y,kernel=.)$predict(xx)}, target=yy)
mbc(Gaussian, Matern32, Matern52,
    evaluator={GauPro_kernel_model$new(X=x,Z=y,kernel=.)$predict(xx)},
    target="yy", inputi={xx <- matrix(runif(nn*d), nn,d);yy <- f(xx)})
mbc(Gaussian, Matern32, Matern52,
    evaluator={GauPro_kernel_model$new(X=x,Z=y,kernel=.)$predict(xx, se=T)},
    target="yy", inputi={xx <- matrix(runif(nn*d), nn,d);yy <- f(xx)}, metric=c('mis90', 'rmse'))
# Use different kinds of model
mbc(KGauss=GauPro_kernel_model$new(X=x,Z=y,kernel=Gaussian),
    Kmat32=GauPro_kernel_model$new(X=x,Z=y,kernel=Matern32),
    Kmat52=GauPro_kernel_model$new(X=x,Z=y,kernel=Matern52),
    GPG=GauPro_Gauss$new(X=x,Z=y), GPGLOO=GauPro_Gauss_LOO$new(X=x,Z=y), times=5,
    inputi={x <- matrix(runif(n*d), n,d);y <- f(x)},
    post=function(x)x$predict(xx,se=T), target=yy, metric=c("rmse","mis90"))

