# Read in data
df <- read.csv("C:\\Users\\colli\\Downloads\\peak_example.csv")

# Transform into reasonable ranges
df$y2 <- df$y/max(df$y)
df$x2 <- (df$x - min(df$x)) / (max(df$x) - min(df$x))

library(GauPro)
# Fit model
gp <- gpkm(df, y2 ~ x2)

# Check fit, make sure it looks reasonable
gp$plot()

# Find x that gives maximum
gpmax <- gp$optimize_fn(fn=function(x) {gp$predict(x)})

# Check the prediction at that point, get an estimate of the std error
gp$predict(gpmax$par$x2, se.fit=T)
