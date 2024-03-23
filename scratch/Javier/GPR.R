library(GauPro)
library(tidyr)
library(fitdistrplus)
# setwd("/Users/javspirata/Documents/Tesis de Gradi/Rstudio")

# datos <- read.table("NMC_ok.csv", header = TRUE, sep = ";", skip= 1)
# datos <- read.table("./scratch/Javier/NMC_ok.csv", header = TRUE, sep = ";", skip= 1)
# datos <- readr::read_csv2(file="./scratch/Javier/NMC_ok.csv")
datos <- readr::read_csv(file="./scratch/Javier/NMC_fixed.csv", skip=1)
# datos$Capacidad <- as.numeric(gsub(",", ".", datos$Capacidad))
plot(datos$ciclo, datos$Capacidad, main = "Capacidad LCO", xlab = "Ciclo", ylab = "Capacidad [%]")
Ciclo <- datos$ciclo
Capacidad <- datos$Capacidad
gp <- GauPro(Ciclo, Capacidad , parallel=FALSE)
plot(Ciclo, Capacidad)
curve(gp$predict(x), add=T, col=2)
curve(gp$predict(x), add=T, col=2)
curve(gp$predict(x)+2*gp$predict(x, se=T)$se, add=T, col=4)
curve(gp$predict(x)-2*gp$predict(x, se=T)$se, add=T, col=4)

if (requireNamespace("MASS", quietly = TRUE)) {
  plot(gp)
}

kern <- Matern52$new(0)
gpk <- GauPro_kernel_model$new(matrix(datos$ciclo, ncol=1), datos$Capacidad , kernel=kern, parallel=FALSE)

gpk$cool1Dplot(
#  n2 = 20,
# nn = 201,
#col2 = "green",
# xlab = "x",
#ylab = "y",
#xmin = NULL,
#xmax = NULL,
# ymin = NULL,
# ymax = NULL,
# gg = TRUE
)
#gpk$plotLOO()

#plot(gpk)

plot(datos$ciclo, datos$Capacidad, main = "LCO Kernel Matern 5/2", xlab = "Ciclo", ylab = "Capacidad [%]")
curve(gpk$predict(x), add=T, col=5)
n <- gp$predict(datos)
plot(gp)
gp$plot1D()



kern <- Matern52$new(0)
gpk <- GauPro_kernel_model$new(matrix(datos$ciclo, ncol=1), datos$Capacidad , kernel=kern, parallel=FALSE)
gpk$predict(datos$ciclo, se.fit = T)

datos %>% tail
gpk$predict(511:1000, se.fit = T)
gpk$plot1D(xmax=1000)

pred <- gpk$predict(1:1000, se.fit=T)
plot(datos$ciclo, datos$Capacidad, main = "LCO Kernel Matern 5/2", xlab = "Ciclo", ylab = "Capacidad [%]", xlim=c(0,1000))
plot(pred[,1], type='l')
plot(datos)

gpklinear <- GauPro_kernel_model$new(trend=trend_LM$new(D=1), matrix(datos$ciclo, ncol=1), datos$Capacidad , kernel=kern, parallel=FALSE)
gpklinear$plot1D(xmax=1000)
