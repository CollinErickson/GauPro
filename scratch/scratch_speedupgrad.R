# Trying to speed grad up
d <- 6
n <- 60
x <- lhs::optimumLHS(n=n,k=d)
f <- TestFunctions::OTL_Circuit
y <- f(x)
gp <- GauPro::GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian)
x1 <- runif(6)
microbenchmark::microbenchmark(gp$grad(x1))
# mean 2838 microsec, median 846 microsec
# after using Kinv down to mean 878, median 670
# after using dC_dx_arma down to mean 166 median 137
pv <- profvis::profvis(replicate(1000, gp$grad(x1)))
pv

gp$kernel$dC_dx(XX = matrix(x1,ncol=d), X=gp$X)
gp$kernel$dC_dx_arma(XX = matrix(x1,ncol=d), X=gp$X)
summary(c(gp$kernel$dC_dx(XX = matrix(x1,ncol=d), X=gp$X)-
          gp$kernel$dC_dx_arma(XX = matrix(x1,ncol=d), X=gp$X)))

microbenchmark::microbenchmark(gp$kernel$dC_dx(XX = matrix(x1,ncol=d), X=gp$X),
                               gp$kernel$dC_dx_arma(XX = matrix(x1,ncol=d), X=gp$X))

pv2 <- profvis::profvis(gp <- GauPro::GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian))
pv2
# Run 10 times
pv2 <- profvis::profvis(replicate(10,GauPro::GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian)))
pv2

# Trying to speed up optimization, solve(C, di) is bottleneck
gp <- GauPro::GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian)
microbenchmark::microbenchmark(gp <- GauPro::GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian), times=5)
# Original:     Mean 2.7088 sec, median 2.714483 sec
# After using Cinv: Mean 1.567 sec, median 1.545 sec
# After removing matmul: Mean 1.27 sec, median 1.26
# After arma gradfuncarray: Mean .683 (1.059), median .684 (.844) on 5 (50) times
