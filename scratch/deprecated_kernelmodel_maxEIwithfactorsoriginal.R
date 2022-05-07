# This was removed from kernel_model.R.
# This ran optim over all factor combinations.
# This can be a lot of optims, so the new one does
# alternating descents on the cts and factor variables.

maxEIwithfactorsorig = function(lower=apply(self$X, 2, min),
                                upper=apply(self$X, 2, max),
                                n0=100, minimize=FALSE, eps=0) {
  stop("dont use this anymore")
  stopifnot(all(lower < upper))
  stopifnot(length(n0)==1, is.numeric(n0), n0>=1)
  # Get factor info
  factorinfo <- find_kernel_factor_dims(self$kernel)
  # Run inner EI over all factor combinations
  stopifnot(length(factorinfo)>0)
  # factordf <- data.frame(index=factorinfo[1]
  factorlist <- list()
  for (i_f in 1:(length(factorinfo)/2)) {
    factorlist[[as.character(factorinfo[i_f*2-1])]] <- 1:factorinfo[i_f*2]
  }
  factordf <- do.call(expand.grid, factorlist)
  # Track best seen in optimizing EI
  bestval <- Inf
  bestpar <- c()
  factorxindex <- factorinfo[(1:(length(factorinfo)/2))*2-1] #factorinfo[[1]]
  for (i_indcomb in 1:prod(factorinfo[(1:(length(factorinfo)/2))*2])) {
    factorxlevel <- unname(unlist(factordf[i_indcomb,])) #i_indcomb #factorinfo[[2]]
    # cat(factorxindex, factorxlevel, "\n")

    # If no non-factor levels, just calculate and compare
    if (length(factorxindex) == self$D) {
      # stop()
      xxinds1 <- c()
      xxinds2 <- c()
      xx <- rep(NA, self$D)
      xx[factorxindex] <- factorxlevel
      optim_out_i_indcomb <- list(par=xx)
      optim_out_i_indcomb$value <- -self$EI(xx, minimize = minimize)
    } else {

      # Otherwise optimize over continuous values
      X0 <- lhs::randomLHS(n=n0, k=self$D)
      X0 <- sweep(X0, 2, upper-lower, "*")
      X0 <- sweep(X0, 2, lower, "+")
      for (j in 1:length(factorxindex)) {
        X0[, factorxindex[j]] <- factorxlevel[j]
      }

      # Calculate EI at these points, use best as starting point for optim
      EI0 <- self$EI(x=X0, minimize=minimize, eps=eps)
      ind <- which.max(EI0)

      # Continuous indexes
      ctsinds <- setdiff(1:self$D, factorxindex)

      # Optimize starting from that point to find input that maximizes EI
      optim_out_i_indcomb <- optim(par=X0[ind, -factorxindex],
                                   lower=lower[-factorxindex], upper=upper[-factorxindex],
                                   # fn=function(xx){ei <- -self$EI(xx); cat(xx, ei, "\n"); ei},
                                   fn=function(xx){
                                     xx2 <- numeric(self$D)
                                     xx2[ctsinds] <- xx
                                     xx2[factorxindex] <- factorxlevel
                                     # xx2 <- c(xx[xxinds1], factorxlevel, xx[xxinds2])
                                     # cat(xx, xx2, "\n")
                                     -self$EI(xx2, minimize = minimize)
                                   },
                                   method="L-BFGS-B")
    }
    if (optim_out_i_indcomb$value < bestval) {
      # cat("new best val", optim_out_i_indcomb$value, bestval, i_indcomb, "\n")
      bestval <- optim_out_i_indcomb$value

      bestpar <- numeric(self$D)
      bestpar[ctsinds] <- optim_out_i_indcomb$par #[xxinds1]
      bestpar[factorxindex] <- factorxlevel
      # bestpar <- c(optim_out_i_indcomb$par[xxinds1], factorxlevel, optim_out_i_indcomb$par[xxinds2])
    }
  }
  stopifnot(length(bestpar) == self$D)
  return(bestpar)
},
