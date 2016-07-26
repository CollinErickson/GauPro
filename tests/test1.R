source("R/corr.R")
gp <- GauPro$new()
gp$fit(matrix(runif(6),3,2),1:3)
gp$pred(matrix(runif(6),3,2))

# 1D test
x <- matrix(seq(0,1,length.out = 9), ncol=1)
y <- sin(2*pi*x)
gp <- GauPro$new()
gp$fit(x,y)
gp$pred(matrix(runif(6),3,2))
curve(gp$pred(x)$me)


# 2D test
n <- 50
x <- matrix(runif(n*2), ncol=2)
f1 <- function(a) {sin(4*pi*a[1]) + sin(6*pi*a[2])}
y <- apply(x,1,f1)
contourfilled::contourfilled.data(x,y)
gp <- GauPro$new()
gp$fit(x,y)
gp$pred(matrix(runif(6),3,2))
contourfilled::contourfilled.func(gp$pred,out.name = 'mean', batchmax = 200)
