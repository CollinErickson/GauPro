
f <- TestFunctions::piston
d <- 7
n <- 60
x <- lhs::randomLHS(n=n,k=d)
y <- f(x)
# y
system.time({gp <- GauPro_kernel_model$new(X=x, Z=y, kernel = Matern52$new(D=d), verbose = 5)})
system.time({gp <- GauPro_kernel_model$new(X=x, Z=y, kernel = Gaussian$new(D=d), verbose = 5)})
plot(gp$pred_LOO(), y)
gp$plotmarginal()
gp$plotmarginal(gp$X[1,])
plot(gp)
gp$plotmarginalrandom()
gp$EI(runif(7))
gp$EI(lhs::randomLHS(n=100, k=ncol(x)))
xmx <- c(1.2770051,  -0.2920814,   0.9825472,  -0.2937785,  -1.3244573,   6.8359251, -11.4165417)
optim(par=xmx, fn=function(xx){ei <- -gp$EI(xx); cat(xx, ei, "\n"); ei})
gp$maxEI()
gp$maxEI(minimize = T)
f(gp$maxEI())
f(gp$maxEI(minimize = T))
gp$maxqEI(5)
f(gp$maxqEI(5))
gp$maxqEI(5, minimize = T)
f(gp$maxqEI(5, minimize = T))

reldiff <- function(a,b) {abs(a-b)/max(abs(c(a,b)))}

