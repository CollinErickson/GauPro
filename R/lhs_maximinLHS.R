lhs_maximinLHS <- function(n, k) {

  if (requireNamespace("lhs", quietly = TRUE)) {
    lhs <- lhs::maximinLHS(n=n, k=k)
  } else {
    message("lhs package not available, using worse option. Please install lhs.")
    # Increasing lhs
    lhs <- (matrix(data=1:n, byrow=F,
                   nrow=n, ncol=k) - 1 +
              matrix(data=runif(n*k),
                     nrow=n, ncol=k)
    ) / n
    # Randomize each column
    for (i in 1:k) {
      lhs[, i] <- lhs[sample(1:n, n, replace=F), i]
    }
  }
  lhs
}
if (F) {
  ceiling(lhs_maximinLHS(10, 3)*10)
}
