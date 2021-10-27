

f <- function(x) {abs(sin(2*pi*x^1.3))^1.3}
# f <- function(x) sin(2*pi*x)
d <- 1
n <- 7
x <- lhs::randomLHS(n=n,k=d)
noisesd <- 1e-16
y <- f(x) + rnorm(n,0, noisesd)
plot(x,y)
# y
system.time({gp <- GauPro_kernel_model$new(X=x, Z=y, kernel = Matern52$new(D=d), verbose = 5)})
# system.time({gp <- GauPro_kernel_model$new(X=x, Z=y, kernel = Gaussian$new(D=d), verbose = 5, restarts = 0)})
plot(gp)
gp$plot1D()
gp$cool1Dplot()
gp$pred(matrix(c(.1,.2,.3,.4,.5), ncol=1), se.fit = T, mean_dist = T)
gp$pred(matrix(c(.1,.2,.3,.4,.5), ncol=1), se.fit = T, mean_dist = F)
curve(gp$EI(matrix(x,ncol=1)) %>% {./(max(.)-min(.))}, add=T, col=3)
gp$EI(matrix(seq(0,1,l=101),ncol=1)) %T>% plot %>% summary
gp$maxqEI(5)

# Run EI
for (i in 1:15) {
  x.ei <- gp$maxEI(lower=0, upper=1, minimize = F)
  gp$update(Xnew=x.ei, Znew=f(x.ei) + rnorm(1,0, noisesd))
  gp$cool1Dplot()
}
