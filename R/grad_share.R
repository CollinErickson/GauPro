#grad_share <- function(fn_gr) {
#  function(x) {
#    out <- fn_gr(x)
#    grad_store <<- out[[2]]
#    out[[1]]
#  }
#}

grad_share <- function(fn_gr) {
  env <- new.env()
  env$fn <- function(x) {
    out <- fn_gr(x)
    env$x_last <- x
    env$fn_val <- out[[1]]
    env$gr_val <- out[[2]]
    #grad_store <<- out[[2]]
    #out[[1]]
    env$fn_val
  }
  env$gr <- function(x = NULL) {
    # Can check if evaluated at same value, but will only slow it down
    #if (!is.null(x) && !any(is.nan(x)) && x != env$x_last) {warning("gr called at different x than fn")}
    env$gr_val
  }
  env
}


fngr <- function(fn_gr, check_all=FALSE, recalculate_indices = 1) {
  env <- new.env()
  env$f <- function(i, check=check_all, recalculate = any(i==recalculate_indices)) {
    function(x=NULL, check_now=check, recalculate_now=recalculate) {
      if (recalculate_now) {
        out <- fn_gr(x)
        env$x_last <- x
        env$out <- out
        out[[1]]
      } else {
        # Can check if evaluated at same value, but will only slow it down
        if (check_now) {
          if (!is.null(x) && !any(is.nan(x)) && x != env$x_last) {
            warning("gr called at different x than fn")
          }
        }
      }
      env$out[[i]]
    }
  }
  env
}

quad_share <- function(x){list(sum(x^4), 4*x^3)}

#grad_share(quad_share)

#ptim_share <- function(par, fngr, ...) {
#  fn <- grad_share(fngr)
#  optim(par=par, fn=fn, gr=function(xx) {grad_store}, ...)
#}
optim_share <- function(par, fngr, ...) {
  env <- grad_share(fngr)
  optim(par=par, fn=env$fn, gr=env$gr, ...)
}
#optim_share(par=c(3, -5), quad_share, method="BFGS")

#lbfgs_share <- function(fngr, vars, method=NULL,...) {
#  fn <- grad_share(fngr)
#  lbfgs::lbfgs(call_eval=fn, call_grad=function(xx) {grad_store}, vars=vars, ...)
#}

lbfgs_share <- function(fngr, vars, method=NULL,...) {
  env <- grad_share(fngr)
  lbfgs::lbfgs(call_eval=env$fn, call_grad=env$gr, vars=vars, ...)
}
#lbfgs_share(vars=c(3, -5), fngr=quad_share)
#parallel::mclapply(1:10, function(i)lbfgs_share(vars=runif(2,-10,10), fngr=quad_share, invisible=1), mc.cores=1)

make_share <- function(func, arg_fn, arg_gr) {
  function(fngr, ...) {
    env <- grad_share(fngr)
    args_list <- list(env$fn, env$gr, ...)
    names(args_list)[1] <- arg_fn
    names(args_list)[2] <- arg_gr
    do.call(what=func, args=args_list)
  }
}
# make_share(lbfgs::lbfgs, 'call_eval', 'call_grad')
# make_share(lbfgs::lbfgs, 'call_eval', 'call_grad')(quad_share, vars=c(5,-4))
