EGO <- R6::R6Class(
  classname = "EGO",
  public = list(
    gp = NULL,
    initialize = function(X, Z, n) {
      self$gp <- GauPro_kernel_model$new(X=X, Z=Z, kernel=kernel)
      self$n <- n
      self$run()
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
