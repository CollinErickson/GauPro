n <- 50
x <- lhs::maximinLHS(n=n, k=2)
y <- TestFunctions::banana(x)
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian)
xx <- matrix(runif(2e5),ncol=2)
gp$pred(XX = xx, se=T)
pv <- profvis::profvis(gp$pred(XX = xx, se=T))
pv


system.time(replicate(10,gp$pred(XX=xx, se=T))) # 23.01 0.06 25.42 before changing df to cbind
system.time(replicate(10,gp$pred(XX=xx, se=T))) # 12.50 0.02 13.12 after changing df to list

microbenchmark::microbenchmark(as.data.frame(tc), as.data.frame.matrix(tc))
microbenchmark::microbenchmark(data.frame(c1,c2,c3), {(list(c1,c2,c3))}, cbind(c1,c2,c3))
