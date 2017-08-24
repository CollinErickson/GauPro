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
if (requireNamespace("MASS", quietly = TRUE)) {
  plot(gp)
}

## ------------------------------------------------------------------------
kern <- Matern52$new(0)
gpk <- GauPro_kernel_model$new(matrix(x, ncol=1), y, kernel=kern, parallel=FALSE)
if (requireNamespace("MASS", quietly = TRUE)) {
  plot(gpk)
}

## ------------------------------------------------------------------------
kern.exp <- Exponential$new(0)
gpk.exp <- GauPro_kernel_model$new(matrix(x, ncol=1), y, kernel=kern.exp, parallel=FALSE)
if (requireNamespace("MASS", quietly = TRUE)) {
  plot(gpk.exp)
}

## ------------------------------------------------------------------------
kern.exp <- Exponential$new(0)
trend.0 <- trend_0$new()
gpk.exp <- GauPro_kernel_model$new(matrix(x, ncol=1), y, kernel=kern.exp, trend=trend.0, parallel=FALSE)
if (requireNamespace("MASS", quietly = TRUE)) {
  plot(gpk.exp)
}

