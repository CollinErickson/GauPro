source("R/corr.R")
gp <- GauPro$new()
gp$fit(matrix(runif(6),3,2),1:3)
gp$pred(matrix(runif(6),3,2))
