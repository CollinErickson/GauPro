corr_gauss <- function(a, b, theta) {#browser()
  exp(-sum(theta * (a-b)^2))
}

corr_gauss_matrix <- function(x, x2=x, theta) {#browser()
  #outer(x,x2, gauss_cor)
  outer(1:nrow(x),1:nrow(x2), Vectorize(function(i,j) corr_gauss(x[i,], x2[j,], theta=theta)))
}
