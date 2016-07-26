library(R6)
GauPro <- R6Class(classname = "GauPro",
  public = list(
    X = NA,
    Z = NA,
    N = NA,
    D = NA,
    corr = "Gauss",
    nug = 1e-6,
    theta = 1,
    fit = function(X, Z) {
      self$X <- X
      self$Z <- Z
      self$N <- nrow(X)
      self$D <- ncol(X)
    },
    pred = function(XX) {#browser()
      if (!is.matrix(XX)) {
        if (self$D == 1) XX <- matrix(XX, ncol=1)
        else if (length(XX) == self$D) XX <- matrix(XX, nrow=1)
        else stop('Predict input should be matrix')
      }
      #covmat <- gauss_cor(c(x, xx))
      kx <- gauss_cor_mat(self$X) + diag(self$nug, self$N)
      kxx <- gauss_cor_mat(XX)
      kx.xx <- gauss_cor_mat(XX, self$X)

      mn = kx.xx %*% solve(kx, self$Z)
      se <- 1 * diag(kxx - kx.xx %*% solve(kx, t(kx.xx)))
      list(mean=mn, se=se)
    }
  ),
  private = list(

  )
)
