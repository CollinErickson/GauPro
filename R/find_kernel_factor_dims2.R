# Change this so Factor/Latent are one group and Ordered are other
# (xindex, nlevels, type)
find_kernel_factor_dims2 <- function (kern) {
  if (("GauPro_kernel_product" %in% class(kern)) || ("GauPro_kernel_sum" %in% class(kern))) {
    return(c(find_kernel_factor_dims2(kern$k1),
             find_kernel_factor_dims2(kern$k2)))
  }
  if (("GauPro_kernel_FactorKernel" %in% class(kern)) ||
      ("GauPro_kernel_LatentFactorKernel" %in% class(kern))) {
    return((c(kern$xindex, kern$nlevels, 0)))
  }
  if (("GauPro_kernel_OrderedFactorKernel" %in% class(kern))) {
    return((c(kern$xindex, kern$nlevels, 1)))
  }
  if (("GauPro_kernel_IgnoreInds" %in% class(kern))) {
    t1 <- find_kernel_factor_dims2(kern$kernel)
    if (is.null(t1)) {
      return(NULL)
    }
    for (i in 1:(length(t1)/2)) {
      # t1[2*i-1] <- t1[2*i-1] + sum(t1[2*i-1] <= kern$ignoreinds)
      t1[2*i-1] <- setdiff((1:(t1[1]+max(kern$ignoreinds))),
                           kern$ignoreinds)[t1[2*i-1]]
    }
    return(t1)
  }
  return(NULL)
}
if (F) {
  k1 <- Gaussian$new(D=2)
  find_kernel_factor_dims(k1)
  find_kernel_factor_dims2(k1)
  k1 <- OrderedFactorKernel$new(D=2, xindex = 2, nlevels = 3)
  find_kernel_factor_dims(k1)
  find_kernel_factor_dims2(k1)
  k1 <- OrderedFactorKernel$new(D=2, xindex = 2, nlevels = 3) *
    FactorKernel$new(D=2, xindex = 1, nlevels = 3)
  find_kernel_factor_dims(k1)
  find_kernel_factor_dims2(k1)

}
