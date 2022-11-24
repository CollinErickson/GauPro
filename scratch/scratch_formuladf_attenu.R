library(dplyr)
x <- datasets::attenu[,c(1,2,4)]
y <- datasets::attenu[,5]
str(x)
str(y)

gp <- GauPro::GauPro_kernel_model$new(X=x, Z=y, kernel='m52')
gp$plot()
gp$plotmarginalrandom()
gp$plotLOO()

gpdf <- GauPro::GauPro_kernel_model$new(datasets::attenu, accel ~ event + mag + dist)
gpdf
gpdf$plotmarginal()
gpdf$plotLOO()
summary(gpdf)
gpdf$importance()
