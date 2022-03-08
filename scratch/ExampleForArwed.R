# Install from my GitHub for most up to date
# devtools::install_github("CollinErickson/GauPro")

# Make fake data
n <- 20
x <- runif(n)
y <- 1.4*x^1.2 - 2.8*sin(2*x) + rnorm(n, 0, 1)
plot(x,y)


# Set parameters
lambda <- 1/4
tau_sq <- 4
nug <- 1

# Load library
library(GauPro)
# Mean/trend function is 0
gp_mean <- GauPro::trend_0$new()
# Kernel with fixed parameters
gp_kernel <- GauPro::Gaussian$new(D=1,
                                  beta=log(1/lambda, 10), beta_est=F,
                                  s2=tau_sq, s2_est=F)
# Fit GP. If nugget is set, there's nothing to estimate
gp <- GauPro::GauPro_kernel_model$new(trend=gp_mean,
                                      kernel=gp_kernel,
                                      X=matrix(x, ncol=1),
                                      Z=y,
                                      nug=nug, nug.est=F)
gp$nug

# Plot what samples look like
gp$cool1Dplot()

# Get samples at specific points
# This does 5 samples at 100 points from 0 to 1.
gp$sample(XX=matrix(seq(0, 1, l=100), ncol=1), 5)
