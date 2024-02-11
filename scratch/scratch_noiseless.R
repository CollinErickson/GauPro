d <- 1
n <- 20
X <- matrix(runif(d*n), ncol=d)
Y <- X[,1]^1.3 - cos(5*sqrt(.2+X[,1]))
cbind(X, Y)
plot(X, Y)

gp <- GauPro_kernel_model$new(X, Y, kernel=Matern52$new(0), nug.max=0, nug.min=0)
gp$plot1D()
gp$nug
cbind(Y, gp$pred(X))
