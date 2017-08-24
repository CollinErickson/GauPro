#' EGO R6 class
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @useDynLib GauPro, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats optim
#' @keywords optimization, Bayesian optimization
#' @return Object of \code{\link{R6Class}} with methods for running EGO.
#' @format \code{\link{R6Class}} object.
#' @examples
#' e1 <- EGO$new(func=sin, n0=10, n=10, d=1)
EGO <- R6::R6Class(
  classname = "EGO",
  public = list(
    gp = NULL,
    func = NULL,
    n0 = NULL,
    n = NULL,
    X = NULL,
    Z = NULL,
    initialize = function(func, n0, n, d) {
      self$func <- func
      self$n0 <- n0
      self$n <- n
      self$d <- d

      self$initial_run()
      self$run()
    },
    initial_run = function() {
      self$X <- lhs::randomLHS(n=n0, k=d)
      self$Z <- apply(self$X, 1, self$func)
      self$gp <- GauPro_kernel_model$new(X=X, Z=Z, kernel=kernel)
    },
    run = function() {
      for (i in 1:n) {
        self$run1()
      }
    },
    run1 = function() {
      # Optimize EI with many start points

      # Select best
      Xbest
      Zbest <- self$func(Xbest)

      # Add to X and Z
      self$X <- rbind(self$X, Xbest)
      self$Z <- c(self$Z, Zbest)

      # Update model
      self$gp$update(Xall=self$X, Zall=self$Z)
    }
  )
)
