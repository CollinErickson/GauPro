EGO <- R6::R6Class(
  classname = "EGO",
  public = list(
    gp = NULL,
    func = NULL,
    n0 = NULL,
    n = NULL,
    initialize = function(func, X, Z, n0, n) {
      self$gp <- GauPro_kernel_model$new(X=X, Z=Z, kernel=kernel)
      self$func <- func
      self$n0 <- n0
      self$n <- n

      self$initial_run()
      self$run()
    },
    initial_run = function() {
      X <- lhs::randomLHS(n=n0, k=)
    },
    run = function() {
      for (i in 1:n) {
        self$run1()
      }
    },
    run1 = function() {

    }
  )
)
