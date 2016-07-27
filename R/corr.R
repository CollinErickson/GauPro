gauss_cor <- function(a, b, theta) {#browser()
  exp(-sum(theta * (a-b)^2))
}

gauss_cor_mat <- function(x, x2=x, theta) {#browser()
  #outer(x,x2, gauss_cor)
  outer(1:nrow(x),1:nrow(x2), Vectorize(function(i,j) gauss_cor(x[i,], x2[j,], theta=theta)))
}
