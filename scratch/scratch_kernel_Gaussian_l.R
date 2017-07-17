# Check numerically that gradient is correct for 1D
# For Gaussian_l kernel
set.seed(0)
n <- 20
x <- matrix(seq(0,1,length.out = n), ncol=1)
f <- Vectorize(function(x) {sin(2*pi*x) + .5*sin(4*pi*x) +rnorm(1,0,.3)})
y <- f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian_l$new(1), parallel=FALSE, verbose=10, nug.est=T)
gp$cool1Dplot()
numDeriv::grad(func = function(x) gp$deviance(params=x[1:2], nuglog=x[3]), x=c(2,1, -6))
gp$deviance_grad(params = c(2,1), nug.update=T, nuglog=-6)
dgc <- function() {
  y <- c(rnorm(2,0,1), rnorm(1,-4,1))
  v1 <- numDeriv::grad(func = function(x) gp$deviance(params=x[1:2], nuglog=x[3]), x=y)
  v2 <- gp$deviance_grad(params = y[1:2], nug.update=T, nuglog=y[3])
  cbind(y, v1, v2, v1-v2, v1/v2)
}
numDeriv::grad(func = function(x) gp$deviance(params=x[1:2], nuglog=x[3]),
               x=c(log(gp$kernel$l,10),gp$kernel$logs2, log(gp$nug,10)))
gp$deviance_grad(params = c(log(gp$kernel$l,10),gp$kernel$logs2), nug.update=T, nuglog=log(gp$nug,10))

# Check dC_dlogl
l <- gp$kernel$l
s2 <- gp$kernel$s2
nug <- 1e-4; lognug=log(nug,10)
eps <- 1e-8
m1 <- (gp$kernel$k(gp$X, l=l+eps, s2=s2) - gp$kernel$k(gp$X, l=l-eps, s2=s2)) / eps / 2
C_nonug <- gp$kernel$k(gp$X, l=l, s2=s2)
C <- C_nonug + diag(nug, nrow(C_nonug))
m2 <- gp$kernel$dC_dparams(params = c(log(l,10), s2), X = gp$X, C = C, C_nonug = C)[[1]][[1]]
c(m1 * l * log(10) -m2) %>% summary
# Check dC_dlogs2
m1 <- (gp$kernel$k(gp$X, l=l, s2=s2+eps) - gp$kernel$k(gp$X, l=l, s2=s2-eps)) / eps / 2
C_nonug <- gp$kernel$k(gp$X, l=l, s2=s2)
C <- C_nonug + diag(nug, nrow(C_nonug))
m2 <- gp$kernel$dC_dparams(params = log(c(l, s2),10), X = gp$X, C = C, C_nonug = C)[[1]][[2]]
c(m1 * s2 * log(10) -m2) %>% summary

numDeriv::grad(func = function(x) gp$deviance(params=x[1:2], nuglog=x[3]),
               x=c(log(c(l,s2),10),log(nug,10)))
gp$deviance_grad(params = log(c(l, s2),10), nug.update=T, nuglog=lognug)



optim(par = c(gp$kernel$logl, gp$kernel$logs2, log(gp$nug,10)),
      fn = function(x) gp$deviance(params = x[1:2], nuglog=x[3]),
      gr = function(x) gp$deviance_grad(params = x[1:2], nuglog=x[3])
        )


# Check 2D
set.seed(0)
n <- 30
x <- lhs::maximinLHS(n=n, k=2)
f <- function(x) {sin(2*pi*x[1]) + .5*sin(4*pi*x[1]) +rnorm(1,0,.03) + x[2]^2}
y <- apply(x, 1, f) #f(x) #sin(2*pi*x) #+ rnorm(n,0,1e-1)
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian_l$new(c(1, 1)), parallel=FALSE, verbose=10, nug.est=T)
ContourFunctions::cf(gp$predict, pts=x)
ContourFunctions::cf(f, pts=x)
numDeriv::grad(func = function(x)gp$deviance(params = x[1:3], nuglog=x[4]), x=c(1,1, 1, -4))
gp$deviance_grad(params = c(1,1,1), nug.update=T, nuglog=-4)
