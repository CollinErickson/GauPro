#' Gaussian process regression model
#'
#' @description
#' Fits a Gaussian process regression model to data.
#'
#' An R6 object is returned with many methods.
#'
#' `gpkm()` is an alias for `GauPro_kernel_model$new()`.
#' For full documentation, see documentation for `GauPro_kernel_model`.
#'
#' Standard methods that work include `plot()`, `summary()`, and `predict()`.
#'
#' @details
#' The default kernel is a Matern 5/2 kernel, but factor/character inputs
#' will be given factor kernels.
#'
#' @export
#' @param X Matrix whose rows are the input points
#' @param Z Output points corresponding to X
#' @param kernel The kernel to use. E.g., Gaussian$new().
#' @param trend Trend to use. E.g., trend_constant$new().
#' @param verbose Amount of stuff to print. 0 is little, 2 is a lot.
#' @param useC Should C code be used when possible? Should be faster.
#' @param useGrad Should the gradient be used?
#' @param parallel Should code be run in parallel? Make optimization
#' faster but uses more computer resources.
#' @param parallel_cores When using parallel, how many cores should
#' be used?
#' @param nug Value for the nugget. The starting value if estimating it.
#' @param nug.min Minimum allowable value for the nugget.
#' @param nug.max Maximum allowable value for the nugget.
#' @param nug.est Should the nugget be estimated?
#' @param param.est Should the kernel parameters be estimated?
#' @param restarts How many optimization restarts should be used when
#' estimating parameters?
#' @param normalize Should the data be normalized?
#' @param optimizer What algorithm should be used to optimize the
#' parameters.
#' @param track_optim Should it track the parameters evaluated
#' while optimizing?
#' @param formula Formula for the data if giving in a data frame.
#' @param data Data frame of data. Use in conjunction with formula.
#' @param ... Not used
gpkm <- function(X, Z,
                 kernel, trend,
                 verbose=0, useC=TRUE, useGrad=TRUE,
                 parallel=FALSE, parallel_cores="detect",
                 nug=1e-6, nug.min=1e-8, nug.max=1e2, nug.est=TRUE,
                 param.est = TRUE, restarts = 0,
                 normalize = FALSE, optimizer="L-BFGS-B",
                 track_optim=FALSE,
                 formula, data,
                 ...) {
  GauPro_kernel_model$new(
    X=X, Z=Z,
    kernel=kernel, trend=trend,
    verbose=verbose, useC=useC, useGrad=useGrad,
    parallel=parallel, parallel_cores=parallel_cores,
    nug=nug, nug.min=nug.min, nug.max=nug.max, nug.est=nug.est,
    param.est = param.est, restarts = restarts,
    normalize = normalize, optimizer=optimizer,
    track_optim=track_optim,
    formula=formula, data=data,
    ...
  )
}
