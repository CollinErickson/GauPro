gauss_cor <- function(a, b, th=.05) {#browser()
  exp(-sum((a-b)^2 / th))
}

gauss_cor_mat <- function(x, x2=x, th=1) {#browser()
  #outer(x,x2, gauss_cor)
  outer(1:nrow(x),1:nrow(x2), Vectorize(function(i,j) gauss_cor(x[i,], x2[j,], th=th)))
}
