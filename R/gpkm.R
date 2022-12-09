#' @rdname GauPro_kernel_model
#' @export
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
