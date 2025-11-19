aq <- airquality
# Add column for day number
aq$daynum <- 1:nrow(aq)
aq <- aq[, -c(5, 6)]
aq
pairs(aq)
# Remove all NA values
aq <- aq[!(is.na(aq$Ozone) | is.na(aq$Solar.R)), ]
pairs(aq)
# Fit Ozone as a factor of the other inputs
gpaq <- gpkm(Ozone ~ ., data=aq, kernel=Matern32, restarts = 20)
gpaq
summary(gpaq)
gpaq$plotmarginalrandom()
gpaq$plotLOO()
gpaq$kernel$k(gpaq$X)
gpaq$kernel$k(gpaq$X, gpaq$X)
