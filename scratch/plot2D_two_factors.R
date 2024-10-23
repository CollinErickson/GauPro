library(dplyr)
library(GauPro)

# numeric x1/x2
n <- 30
x1 <- runif(n)
x2 <- runif(n)
y <- x1*x2 + .02*rnorm(length(x1))
df <- data.frame(x1, x2, y)
df
pairs(df)
gp <- gpkm(y ~ x1 + x2, df)
gp
gp$plot2D()


# factor x1/x2
x1 <- sample(c('a', 'b'), 10, T)
x2 <- sample(c('c', 'd'), 10, T)
y <- as.numeric((x1 == 'a') & (x2 == 'c')) + rnorm(length(x1))
df <- data.frame(x1, x2, y)
df
gp <- gpkm(y ~ x1 + x2, df,
           kernel=k_GowerFactorKernel(D=1, nlevels=2, xindex=1)*
             k_GowerFactorKernel(D=2, nlevels=2, xindex=2))
gp
gp$plot2D()
