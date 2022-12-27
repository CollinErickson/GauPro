
find_kernel_cts_dims <- function (kern) {
  if (("GauPro_kernel_product" %in% class(kern)) || ("GauPro_kernel_sum" %in% class(kern))) {
    return(sort(unique(c(find_kernel_cts_dims(kern$k1),
                         find_kernel_cts_dims(kern$k2)))))
  }
  if (("GauPro_kernel_FactorKernel" %in% class(kern)) ||
      ("GauPro_kernel_LatentFactorKernel" %in% class(kern)) ||
      ("GauPro_kernel_OrderedFactorKernel" %in% class(kern)) ||
      ("GauPro_kernel_GowerFactorKernel" %in% class(kern))) {
    return(NULL)
  }
  if (("GauPro_kernel_IgnoreInds" %in% class(kern))) {
    t1 <- find_kernel_cts_dims(kern$kernel)
    if (is.null(t1)) {
      return(NULL)
    }
    # for (i in 1:length(t1)) {
    #   # t1[2*i-1] <- t1[2*i-1] + sum(t1[2*i-1] <= kern$ignoreinds)
    #   # t1[2*i-1] <- setdiff((1:(t1[1]+max(kern$ignoreinds))),
    #   #                      kern$ignoreinds)[t1[2*i-1]]
    #   t1[i] <- t1[i] + sum(??? < t1[i])
    # }
    a <- 1:(max(kern$ignoreinds, length(kern$ignoreinds) + length(t1)))
    b <- a[-kern$ignoreinds]
    d <- b[t1]
    return(d)
  }
  if ("GauPro_kernel_White" %in% class(kern)) {
    return(NULL)
  }
  if (!("GauPro_kernel" %in% class(kern))) {
    stop("kern isn't a GauPro_kernel")
  }
  # All other kernels are continuous over all dimensions
  return(1:kern$D)
}
if (F) {
  k1 <- Gaussian$new(D=2)
  find_kernel_cts_dims(k1)
  k1 <- OrderedFactorKernel$new(D=2, xindex = 2, nlevels = 3)
  find_kernel_cts_dims(k1)
  k3 <- IgnoreIndsKernel$new(ignoreinds = 3:4, Gaussian$new(D=2)) *
    LatentFactorKernel$new(D=4, nlevels = 2, latentdim = 1, xindex = 3) *
    LatentFactorKernel$new(D=4, nlevels = 4, latentdim = 2, xindex = 4)
  find_kernel_cts_dims(k3)
  k4 <- IgnoreIndsKernel$new(ignoreinds = 1:2, Gaussian$new(D=2))
  find_kernel_cts_dims(k4)
  k4 <- IgnoreIndsKernel$new(ignoreinds = c(1,3), Gaussian$new(D=2))
  find_kernel_cts_dims(k4)
}
