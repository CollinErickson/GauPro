
n <- 12
x <- matrix(seq(0,1,length.out = n), ncol=1)
f <- function(x) {sin(2*pi*x)}
y <- f(x) + rnorm(n,0,1e-2)

gp <- GauPro_Gauss_LOO$new(X=x, Z=y, parallel=FALSE)

curve(f, ylim=c(-2,2))
curve(gp$predict(x), add=T, col=2)
curve(gp$predict(x) + 2 * gp$predict(x,se=T)$se, add=T, col=3)
curve(gp$predict(x) - 2 * gp$predict(x,se=T)$se, add=T, col=3)
points(x, y)

curve(gp$tmod$predict(x))
points(gp$X, gp$tmod$Z)

# Try function with big errors at middle
n <- 13
x <- runif(n)
f <- function(x) {1 / sqrt(1 + exp(-25*(x-.5)))}
y <- f(x) + rnorm(n,0,1e-3)

gp <- GauPro_Gauss_LOO$new(X=x, Z=y, parallel=FALSE)

curve(f, ylim=c(-.05,1.1))
curve(gp$predict(x), add=T, col=2)
curve(gp$predict(x) + 2 * gp$predict(x,se=T)$se, add=T, col=3)
curve(gp$predict(x) - 2 * gp$predict(x,se=T)$se, add=T, col=3)
points(x, y)


# Try function with big errors at middle
n <- 20
x <- runif(n)
f <- function(x) {(2*x) %% 1}
y <- f(x) + rnorm(n,0,1e-1)

gp <- GauPro_Gauss_LOO$new(X=x, Z=y, parallel=FALSE)

curve(f, ylim=c(-.05,1.1))
curve(gp$predict(x), add=T, col=2)
curve(gp$predict(x) + 2 * gp$predict(x,se=T)$se, add=T, col=3)
curve(gp$predict(x) - 2 * gp$predict(x,se=T)$se, add=T, col=3)
points(x, y)
curve(gp$tmod$predict(x))


# Try function that has bigger errors to right
n <- 20
x <- runif(n)
f <- function(x) {x}
y <- f(x) + rnorm(n,0,x/1)

gp <- GauPro_Gauss_LOO$new(X=x, Z=y, parallel=FALSE)

curve(f, ylim=c(-.05,1.1))
curve(gp$predict(x), add=T, col=2)
curve(gp$predict(x) + 2 * gp$predict(x,se=T)$se, add=T, col=3)
curve(gp$predict(x) - 2 * gp$predict(x,se=T)$se, add=T, col=3)
points(x, y)
curve(gp$tmod$predict(x))
points(x, gp$tmod$Z)
cbind(x, y, gp$pred_LOO(se=T))

# Check banana function
n <- 50
x <- lhs::maximinLHS(n=n, k=2)
y <- TestFunctions::banana(x)
gp <- GauPro_Gauss_LOO$new(X=x, Z=y, parallel=FALSE)
ContourFunctions::cf(gp$predict, batchmax=Inf, pts=x)
ContourFunctions::cf(function(x) gp$predict(x, se.fit=T)$se, batchmax=Inf, pts=x)
ContourFunctions::cf(gp$tmod$predict, batchmax=Inf, pts=x)

# Turn LOO off to see difference
gp$use_LOO <- FALSE
ContourFunctions::cf(gp$predict, batchmax=Inf, pts=x)
ContourFunctions::cf(function(x) gp$predict(x, se.fit=T)$se, batchmax=Inf, pts=x)
