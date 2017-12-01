# Testing update_fast
n <- 50
x <- lhs::maximinLHS(n=n, k=2)
y <- TestFunctions::banana(x)
x2 <- matrix(c(.5,.5), 1,2)
y2 <- TestFunctions::banana(x2)
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian)
gp$update_fast(Xnew = x2, Znew = y2)

microbenchmark::microbenchmark({
  add20=for (i in 1:20) {
    x3 <- matrix(runif(2), 1, 2)
    y3 <- TestFunctions::banana(x3)
    gp$update_fast(Xnew = x3, Znew = y3)
    # gp$update(Xnew = x3, Znew = y3)
    # gp$update(Xnew = x3, Znew = y3, no_update = T)
  }}, times=1
)
