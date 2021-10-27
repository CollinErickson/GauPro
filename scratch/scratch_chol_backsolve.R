n <- 500
d <- 1
kern <- GauPro::Gaussian$new(D=d, beta=2)
X <- matrix(runif(n*d), ncol=d)
kmat <- kern$k(x=X)
diag(kmat) <- diag(kmat) + 1e-4
kchol <- chol(kmat)
t(kchol) %*% kchol # equals kmat
kinv <- chol2inv(kchol)
Z <- rnorm(n)
kinv %*% Z
# backsolve(kchol, backsolve(kchol, Z), transpose = T)
# forwardsolve(kchol, backsolve(kchol, Z), transpose = T)
# This is it
backsolve(kchol, backsolve(kchol, Z, transpose = T))
plot(kinv %*% Z, backsolve(kchol, backsolve(kchol, Z, transpose = T)))
# backsolves are way faster than chol2inv
microbenchmark::microbenchmark(
  chol2inv(kchol),
  backsolve(kchol, backsolve(kchol, Z, transpose = T))
)
# Can cut time in half by using backsolves
microbenchmark::microbenchmark(
  {kchol <- chol(kmat); chol2inv(kchol)},
  {kchol <- chol(kmat); backsolve(kchol, backsolve(kchol, Z, transpose = T))}, times=100
)
