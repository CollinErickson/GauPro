EI_tdist <- function(target, y, s, df, minimize=TRUE) {
  if (!minimize) {
    target <- -target
    y <- -y
  }
  z <- (target - y) / s
  (target - y) * pt(z, df) + df / (df - 1) * (1 + z^2/df) * s * dt(z, df)
}
if (F) {
  EI_tdist(0, 1, 1, 10, F)
  EI_tdist(0, 1, 1, 10, T)
  # Vary df
  curve(EI_tdist(0, 1, 1, x, T), 3, 30)
  curve(EI_tdist(0, 1, 1, x, F), 3, 30)
  # Vary mean
  curve(EI_tdist(0, x, 1, 10, T), -5, 5)
  curve(EI_tdist(0, x, 1, 10, F), -5, 5)
  # Vary s
  curve(EI_tdist(0, 1, x, 10, T), .01, 3)
  curve(EI_tdist(0, 1, x, 10, T), .01, 3)
}
