convert_X_with_formula_back <- function(gpdf, x) {

  par <- x

  if (is.matrix(par)) {
    pardf <- as.data.frame(par)
  } else if (is.numeric(par)) {
    pardf <- as.data.frame(matrix(par, nrow=1))
  }
  colnames(pardf) <- colnames(gpdf$X)
  pardf

  # Convert factor indexes back to factor
  for (i in seq_along(gpdf$convert_formula_data$factors)) {
    pardf[[gpdf$convert_formula_data$factors[[i]]$index]] <-
      gpdf$convert_formula_data$factors[[i]]$levels[
        pardf[[gpdf$convert_formula_data$factors[[i]]$index]]]
  }
  # Convert char indexes back to char
  for (i in seq_along(gpdf$convert_formula_data$chars)) {
    pardf[[gpdf$convert_formula_data$chars[[i]]$index]] <-
      gpdf$convert_formula_data$chars[[i]]$vals[
        pardf[[gpdf$convert_formula_data$chars[[i]]$index]]]
  }

  pardf
}
