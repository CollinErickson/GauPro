## ----echo=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(fig.width=7, fig.height=4) 

set.seed(0)

## ----libraryGauPro------------------------------------------------------------
library(GauPro)

## ----fitsine------------------------------------------------------------------
n <- 12
x <- seq(0, 1, length.out = n)
y <- sin(6*x^.8) + rnorm(n,0,1e-1)
gp <- gpkm(x, y)

## ----plotsine-----------------------------------------------------------------
gp$plot1D()

## ----fit_dm-------------------------------------------------------------------
library(ggplot2)
diamonds_subset <- diamonds[sample(1:nrow(diamonds), 60), ]
dm <- gpkm(price ~ carat + cut + color + clarity + depth,
           diamonds_subset)

## ----summary_dm---------------------------------------------------------------
summary(dm)

## ----plot_dm------------------------------------------------------------------
plot(dm)

## ----diamond_construct_kernel-------------------------------------------------
cts_kernel <- k_IgnoreIndsKernel(k=k_PowerExp(D=2), ignoreinds = c(2,3,4))
factor_kernel2 <- k_OrderedFactorKernel(D=5, xindex=2, nlevels=nlevels(diamonds_subset[[2]]))
factor_kernel3 <- k_OrderedFactorKernel(D=5, xindex=3, nlevels=nlevels(diamonds_subset[[3]]))
factor_kernel4 <- k_GowerFactorKernel(D=5, xindex=4, nlevels=nlevels(diamonds_subset[[4]]))

# Multiply them
diamond_kernel <- cts_kernel * factor_kernel2 * factor_kernel3 * factor_kernel4

## ----diamond_construct_kernel_fit---------------------------------------------
dm <- gpkm(price ~ carat + cut + color + clarity + depth,
           diamonds_subset, kernel=diamond_kernel)
dm$plotkernel()

## ----combine seed, include=F--------------------------------------------------
set.seed(99)

## ----combine_periodic---------------------------------------------------------
x <- 1:20
y <- sin(x) + .1*x^1.3
combo_kernel <- k_Periodic(D=1) * k_Matern52(D=1)
gp <- gpkm(x, y, kernel=combo_kernel, nug.min=1e-6)
gp$plot()

## ----oldvignettedata----------------------------------------------------------
x <- seq(0,1,l=10)
y <- abs(sin(2*pi*x))^.8
ggplot(aes(x,y), data=cbind(x,y)) +
  geom_point()

## ----oldvignettedata_plot-----------------------------------------------------
ggplot(aes(x,y), data=cbind(x,y)) +
    geom_point() +
    stat_smooth(method='lm')

## ----oldvignettedata_gpkm-----------------------------------------------------
library(GauPro)
gp <- gpkm(x, y, kernel=k_Gaussian(D=1), parallel=FALSE)

## ----oldvignettedata_plot1D---------------------------------------------------
gp$plot1D()

## ----oldvignettedata_cool1Dplot-----------------------------------------------
if (requireNamespace("MASS", quietly = TRUE)) {
  gp$cool1Dplot()
}

## ----oldvignettedata_maternplot-----------------------------------------------
kern <- k_Matern52(D=1)
gpk <- gpkm(matrix(x, ncol=1), y, kernel=kern, parallel=FALSE)
if (requireNamespace("MASS", quietly = TRUE)) {
  plot(gpk)
}

## ----oldvignettedata_exponentialplot------------------------------------------
kern.exp <- k_Exponential(D=1)
gpk.exp <- gpkm(matrix(x, ncol=1), y, kernel=kern.exp, parallel=FALSE)
if (requireNamespace("MASS", quietly = TRUE)) {
  plot(gpk.exp)
}

## ----oldvignettedata_trendplot------------------------------------------------
kern.exp <- k_Exponential(D=1)
trend.0 <- trend_0$new()
gpk.exp <- gpkm(matrix(x, ncol=1), y, kernel=kern.exp, trend=trend.0, parallel=FALSE)
if (requireNamespace("MASS", quietly = TRUE)) {
  plot(gpk.exp)
}

