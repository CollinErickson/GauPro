find_kernel_factor_dims <- function (kern) {
  if (("GauPro_kernel_product" %in% class(kern)) || ("GauPro_kernel_sum" %in% class(kern))) {
    return(c(find_kernel_factor_dims(kern$k1),
             find_kernel_factor_dims(kern$k2)))
  }
  if (("GauPro_kernel_FactorKernel" %in% class(kern)) || ("GauPro_kernel_OrderedFactorKernel" %in% class(kern))) {
    return((c(kern$xindex, kern$nlevels)))
  }
  return(NULL)
}
