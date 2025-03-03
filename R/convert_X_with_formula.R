#' @importFrom stats model.frame
convert_X_with_formula <- function(X, convert_formula_data, formula) {
  stopifnot(is.data.frame(X))
  data <- X
  # X might not have response var
  stopifnot(length(as.character(formula)) == 3)
  z_name <- as.character(formula)[2]
  if (!(z_name %in% colnames(data))) {
    data[[z_name]] <- 1 # Can't use NA
  }
  modfr <- model.frame(formula = formula, data = data)
  Xdf <- modfr[,2:ncol(modfr)]

  # Convert factor columns to integer
  # for (i in 1:ncol(Xdf)) {
  #   if (is.factor(Xdf[, i])) {
  #     # Check that levels match
  #
  #     convert_formula_data$factors[[
  #       length(convert_formula_data$factors)+1
  #     ]] <- list(index=i,
  #                levels=levels(Xdf[[i]]))
  #     Xdf[[i]] <- as.integer(Xdf[[i]])
  #   }
  # }
  factorinds <- sapply(
    convert_formula_data$factors,
    function(li) {li$index}
  )
  for (iii in seq_along(factorinds)) {
    i <- factorinds[iii]
    # Check that levels match
    # Convert
    if (is.factor(Xdf[[i]])) {
      Xdf[[i]] <- as.integer(Xdf[[i]])
    } else {
      # User can give in character of the level instead of proper factor
      Xdf[[i]] <- sapply(
        Xdf[[i]],
        function(x) {
          which(x == convert_formula_data$factors[[iii]]$levels)
        })
    }
  }
  # # Convert char columns to integer
  # for (i in 1:ncol(Xdf)) {
  #   if (is.character(Xdf[, i])) {
  #     convert_formula_data$chars[[
  #       length(convert_formula_data$chars)+1
  #     ]] <- list(index=i,
  #                vals=sort(unique(Xdf[[i]])))
  #     Xdf[[i]] <- sapply(Xdf[[i]],
  #                        function(x) {
  #                          which(x==convert_formula_data$chars[[
  #                            length(convert_formula_data$chars)
  #                          ]]$vals)
  #                        })
  #   }
  # }

  charinds <- sapply(
    convert_formula_data$chars,
    function(li) {li$index}
  )
  for (ind in seq_along(charinds)) {
    i <- charinds[ind]
    # Check that levels match
    stopifnot(Xdf[[i]] %in% convert_formula_data$chars[[ind]]$vals)
    # Convert
    Xdf[[i]] <- sapply(Xdf[[i]],
                       function(x) {
                         which(x==convert_formula_data$chars[[ind]]$vals)
                       })
  }
  X <- as.matrix(Xdf)
  X
}


if (F) {

  n <- 30
  tdf <- data.frame(a=runif(n), b=runif(n), c=factor(sample(5:6,n,T)), d=runif(n), e=sample(letters[1:3], n,T))
  tdf$z <- with(tdf, a+a*b+b^2)
  tdf[[3]]
  as.integer(tdf[[3]])
  tdf
  gpf <- GauPro_kernel_model$new(X=tdf, Z=z ~ a + b + c + e, kernel='gauss')
  predict(gpf, tdf)
}
