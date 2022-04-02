compM <- comparer::mbc(
  times=10,
  R={
    useCM <<- F
    gp1 <- GauPro_kernel_model$new(
      X=Xmat, Z=y, verbose=3, kernel=Exponential
    )
    gp1$deviance()
  },
  C={
    useCM <<- T
    gp1 <- GauPro_kernel_model$new(
      X=Xmat, Z=y, verbose=3, kernel=Exponential
    )
    gp1$deviance()
  }
)
compM

compMb <- comparer::mbc(
  times=5,
  R={
    useCM <<- F
    gp1 <- GauPro_kernel_model$new(
      X=Xmat, Z=y, verbose=3,
      kernel=IgnoreIndsKernel$new(ignoreinds = 3:4, Exponential$new(D=2)) *
        LatentFactorKernel$new(D=4, nlevels = 2, latentdim = 1, xindex = 3) *
        LatentFactorKernel$new(D=4, nlevels = 4, latentdim = 2, xindex = 4)
    )
    gp1$deviance()
  },
  C={
    useCM <<- T
    gp1 <- GauPro_kernel_model$new(
      X=Xmat, Z=y, verbose=3,
      kernel=IgnoreIndsKernel$new(ignoreinds = 3:4, Exponential$new(D=2)) *
        LatentFactorKernel$new(D=4, nlevels = 2, latentdim = 1, xindex = 3) *
        LatentFactorKernel$new(D=4, nlevels = 4, latentdim = 2, xindex = 4)
    )
    gp1$deviance()
  }
)
compMb
