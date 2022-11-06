# CBE 5/6/2021
# When all optim restarts fail, it sticks with initial values.
# It'd be better if optim could track best value along the way,
# and then if there is an error, return the best value seen so far
# instead of throwing an error.
# In this file optim_share2 improves optim_share to do exactly this.
# There is nearly no increase in runtime on the example I checked.

if (F) {
  optim(par=2, fn=function(x) {x^2}, method="L-BFGS-B")
  optim(par=2, fn=function(x) {x^2}, gr=function(x){2*x}, method="L-BFGS-B")
  optim(par=2, fn=function(x) {ifelse(x>1,x^2, NA)}, gr=function(x){2*x}, method="L-BFGS-B")
  optim_share(par=2, fngr=function(x) {list(x^2, 2*x)}, method="L-BFGS-B")
  optim_share(par=2, fngr=function(x) {list(ifelse(x>1,x^2, NA), 2*x)}, method="L-BFGS-B")
}

optim_share2 <- function(par, fngr, ...) {
  env <- grad_share(fngr)
  bestx <- 0
  besty <- Inf
  # iter <- 0
  neval <- 0
  f1 <- function(x) {
    neval <<- neval + 1
    out <- env$fn(x)
    # iter <<- iter + 1
    # print(c(iter, round(x,4), out, round(bestx,4), besty))
    if (!is.na(out) && out < besty) {
      bestx <<- x
      besty <<- out
    }
    out
  }
  # optim(par=par, fn=env$fn, gr=env$gr, ...)
  optim_out <- try({
    optim(par=par, fn=f1, gr=env$gr, ...)
    # optim(par=par, fn=f1, gr=env$gr, control=list(factr=1e11), ...)
  }, silent = TRUE)
  if (inherits(optim_out, "try-error")) {
    # print('try-error')
    if (is.infinite(besty)) {
      stop(paste0("optim_share2 error on starting params: ",
                  attr(optim_out, 'condition')$message))
    } else {
      return(list(par=bestx,
                  value=besty,
                  counts=c('function'=neval, 'gradient'=NA),
                  convergence=NA,
                  message="FAILED OPTIM: USING BEST SEEN, SEE optim_share2 FOR DETAILS"))
    }
  }
  optim_out
}
if (F) {
  optim_share2(par=2, fngr=function(x) {list(ifelse(x>.1,x^2, NA), 2*x)}, method="L-BFGS-B")
  optim_share2(par=2, fngr=function(x) {list(ifelse(x>.1,x^4, NA), 4*x^3)}, method="L-BFGS-B")
}

if (F) {
  microbenchmark::microbenchmark(os2={use_optim_share2 <- TRUE; k2 <- FactorKernel$new(D=2, nlevels=3, xind=1); gp <- GauPro_kernel_model$new(X=X, Z=Z, kernel = k2, verbose = 5)},
                                 os1={use_optim_share2 <- F; k2 <- FactorKernel$new(D=2, nlevels=3, xind=1); gp <- GauPro_kernel_model$new(X=X, Z=Z, kernel = k2, verbose = 5)}, times=10)

  use_optim_share2 <- TRUE
  k2 <- FactorKernel$new(D=2, nlevels=3, xind=1)
  gp <- GauPro_kernel_model$new(X=X, Z=Z, kernel = k2, verbose = 5)
}
