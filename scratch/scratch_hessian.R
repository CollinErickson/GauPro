
# 2D surface to test Hessian
n <- 40
x <- matrix(runif(n*2), ncol=2)
f1 <- function(a) {sin(2*pi*a[1]) + sin(6*pi*a[2])}
#f1 <- TestFunctions::branin
#f1 <- TestFunctions::RFF_get(D=2)
y <- apply(x,1,f1) + rnorm(n,0,.01)
system.time(cf::cf_data(x,y))
gp <- GauPro(x,y, verbose=2);gp$theta
system.time(cf::cf_func(gp$pred, pts=x))
# They give same numerical answer
gp$hessian(c(.2,.75))
numDeriv::hessian(gp$predict, c(.2,.75))

max_eigen <- function(x) {
  evals <- eigen(gp$hessian(x), symmetric = T, only.values = T)
  maxeval <- evals$val[which.max(abs(evals$val))]
  maxeval
}
cf::cf(max_eigen, batchmax=1, n=40)
