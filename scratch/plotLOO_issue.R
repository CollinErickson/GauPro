# https://github.com/CollinErickson/GauPro/issues/3
# Tried to check his issue, but found a different issue.

n <- 50
f <- function(a,b,c) {a+b*c+100}
smooth_x_rt <- data.frame(a=runif(n), b=rnorm(n), c=rexp(n))
results <- data.frame(RT_2=with(smooth_x_rt, f(a,b,c)), RT=rnorm(n, sd=.1))

gp <- GauPro::gpkm(smooth_x_rt, results$RT_2 - results$RT, kernel = "matern52",
             parallel = TRUE, normalize = TRUE, verbose = 0)
gp$plot()
gp$plotLOO()


n2 <- 1000
x <- data.frame(a=runif(n2), b=rnorm(n2), rexp(n))
predx <- gp$pred(x)
fx <- f(x[,1], x[,2], x[,3])
plot(predx, fx); abline(a=0, b=1)

