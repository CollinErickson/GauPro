# Ordered is way slower than latent

library(dplyr)
n <- 163*2
# Non-ordered data
xdf <- tibble(
  a=rnorm(n),
  b=runif(n),
  c=sample(letters[1:5], n, T),
  d=factor(sample(letters[6:9], n, T)),
  e=rexp(n),
  # f=rnorm(n),
  z=a*b + a^2*ifelse(c %in% c('a', 'b'), 1, .5) +
    e*ifelse(d %in% c('g','h'), 1, -1) +
    ifelse(paste0(d,e) %in% c('af', 'ah', 'cf', 'cg', 'ci'),4,0) +
    rnorm(n, 1e-3)
)
# Ordered data
xdf <- tibble(
  a=rnorm(n),
  b=runif(n),
  c=sample(letters[1:5], n, T),
  d=ordered(sample(letters[6:9], n, T)),
  e=rexp(n),
  # f=rnorm(n),
  z=a*b + a^2*ifelse(c %in% c('a', 'b'), 1, .5) +
    e*ifelse(d %in% c('g','h'), 1, -1) +
    ifelse(paste0(d,e) %in% c('af', 'ah', 'cf', 'cg', 'ci'),4,0) +
    rnorm(n, 1e-3)
)

# Time fit
system.time(expect_error(gpdf <- GauPro_kernel_model$new(z ~ ., data=xdf, kernel='gauss'), NA))
# system.time(expect_error(gpdf <- GauPro_kernel_model$new(z ~ ., data=xdf, kernel='matern52'), NA))

# Latent
# .83, 1.11
# Ordered
# 15.3, 12.7
# Ordered after adding k
# 5.0, 3.0, 4.9
# After adding dc_dparams
# 1.4, 1.8, 1.1, 1.0, 1.5, 1.2
# .7, .75, .83, .77

# pv shows most of time is in k/kone since it does outer
