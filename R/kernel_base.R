# Kernels should implement:
# k kernel function for two vectors
# update_params
# get_optim_functions: return optim.func, optim.grad, optim.fngr
# param_optim_lower - lower bound of params
# param_optim_upper - upper
# param_optim_start - current param values
# param_optim_start0 - some central param values that can be used for optimization restarts
# param_optim_jitter - how to jitter params in optimization

# Suggested
# deviance
# deviance_grad
# deviance_fngr
# grad



#' Kernel R6 class
#'
#' @docType class
#' @importFrom R6 R6Class
# @export
#' @useDynLib GauPro, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats optim
# @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @field D Number of input dimensions of data
#' @examples
#' #k <- GauPro_kernel$new()
GauPro_kernel <- R6::R6Class(classname = "GauPro_kernel",
  public = list(
    D = NULL
    # k_diag = function(x) {
    #   if (is.matrix(x)) {rep(self$s2, nrow(x))}
    #   else if (length(x) == self$D) {self$s2}
    #   else {stop("Error in k_diag #4928")}
    # }
  ),
  private = list(

  )
)
