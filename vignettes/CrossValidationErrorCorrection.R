## ---- results='hide', echo=FALSE-----------------------------------------
set.seed(0)

## ------------------------------------------------------------------------
n <- 200
m1 <- matrix(runif(n*n),ncol=n)
b1 <- runif(n)
if (requireNamespace("microbenchmark", quietly = TRUE)) {
  microbenchmark::microbenchmark(solve(m1, b1), m1 %*% b1)
}

## ------------------------------------------------------------------------
set.seed(0)
corr <- function(x,y) {exp(sum(-30*(x-y)^2))}
n <- 200
d <- 2
X <- matrix(runif(n*d),ncol=2)
R <- outer(1:n,1:n, Vectorize(function(i,j) {corr(X[i,], X[j,])}))
Rinv <- solve(R)
A <- R[-n,-n]
Ainv <- solve(A)
E <- Rinv[-n, -n]
b <- R[n,-n]
g <- Rinv[n,-n]
Ainv_shortcut <- E + E %*% b %*% g / (1-sum(g*b))
summary(c(Ainv - Ainv_shortcut))
if (requireNamespace("microbenchmark", quietly = TRUE)) {
  microbenchmark::microbenchmark(solve(A), E + E %*% b %*% g / (1-sum(g*b)))
}

