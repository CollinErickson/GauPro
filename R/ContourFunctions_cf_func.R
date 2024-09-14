ContourFunctions_cf_func <- function(...) {
  if (requireNamespace("ContourFunctions")) {
    ContourFunctions::cf_func(...)
  } else {
    message(paste0("The R package ContourFunctions is not available. ",
                   "Please install and try again."))
    return(NULL)
  }
}
