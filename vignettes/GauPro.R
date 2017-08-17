## ---- echo=FALSE---------------------------------------------------------
set.seed(0)

## ------------------------------------------------------------------------
x <- seq(0,1,l=10)
y <- abs(sin(2*pi*x))^.8
plot(x, y)

## ------------------------------------------------------------------------
lm_mod <- lm(y ~ x)
plot(x, y)
abline(a=lm_mod$coef[1], b=lm_mod$coef[2], col='red')

## ------------------------------------------------------------------------
library(GauPro)
gp <- GauPro(x, y, parallel=FALSE)

## ------------------------------------------------------------------------
plot(x, y)
curve(gp$predict(x), add=T, col=2)

## ------------------------------------------------------------------------
plot(x, y)
curve(gp$predict(x), add=T, col=2)
curve(gp$predict(x)+2*gp$predict(x, se=T)$se, add=T, col=4)
curve(gp$predict(x)-2*gp$predict(x, se=T)$se, add=T, col=4)

## ------------------------------------------------------------------------
plot(gp)

## ------------------------------------------------------------------------
kern <- Matern52$new(0)
gpk <- GauPro_kernel_model$new(matrix(x, ncol=1), y, kernel=kern)
plot(gpk)

## ------------------------------------------------------------------------
kern <- Exponential$new(0)
gpk <- GauPro_kernel_model$new(matrix(x, ncol=1), y, kernel=kern)
plot(gpk)

