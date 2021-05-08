
f <- TestFunctions::piston
d <- 7
n <- 60
x <- lhs::randomLHS(n=n,k=d)
y <- f(x)
# y
system.time({gp <- GauPro_kernel_model$new(X=x, Z=y, kernel = Matern52$new(D=7), verbose = 5)})
system.time({gp <- GauPro_kernel_model$new(X=x, Z=y, kernel = Gaussian$new(D=7), verbose = 5)})
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

reldiff <- function(a,b) {abs(a-b)/max(abs(c(a,b)))}


# Plot marginal random averages
X <- lhs::randomLHS(n=151, k=ncol(gp$X))
X2 <- sweep(X, 2, apply(gp$X, 2, max) - apply(gp$X, 2, min), "*")
X3 <- sweep(X2, 2, apply(gp$X, 2, min), "+")
X3pred <- gp$pred(X3, se.fit = T)
X3pred$irow <- 1:nrow(X3pred)
head(X3pred)
X4 <- dplyr::inner_join(X3pred, tidyr::pivot_longer(cbind(as.data.frame(X),irow=1:nrow(X)), cols=1:7), "irow")
head(X4)
X4$upper <- X4$mean + 2*X4$se
X4$lower <- X4$mean - 2*X4$se
ggplot(X4, aes(value, mean)) + facet_wrap(.~name) +
  # geom_point(aes(y=upper), color="green") +
  geom_segment(aes(y=upper, yend=lower, xend=value), color="green", size=2) +
  geom_point()

EI = function(self, x, minimize=FALSE, eps=.01) {
  stopifnot(is.vector(x), length(x) == ncol(self$X))
  fxplus <- max(self$Z)
  pred <- self$pred(x, se.fit=T)
  Ztop <- pred$mean - fxplus - eps
  Z <- Ztop / pred$se
  if (pred$se <= 0) {return(0)}
  (Ztop) * pnorm(Z) + pred$se * dnorm(Z)

}
EI(gp, gp$X[2,])
EI(gp, runif(7))
optim(par=runif(7), fn=function(xx){ei <- -predict(gp, xx); cat(xx, ei, "\n"); ei})
optim(par=runif(7), fn=function(xx){ei <- -EI(gp, xx); cat(xx, ei, "\n"); ei})
xmx <- c(1.2770051,  -0.2920814,   0.9825472,  -0.2937785,  -1.3244573,   6.8359251, -11.4165417)
optim(par=xmx, fn=function(xx){ei <- -EI(gp, xx); cat(xx, ei, "\n"); ei})
