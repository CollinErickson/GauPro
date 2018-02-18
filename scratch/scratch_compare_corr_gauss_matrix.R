d <- 6
x1 <- matrix(runif(100*d), ncol=d)
x2 <- matrix(runif(1e4*d), ncol=d)
th <- c(.3,3.3)
t1 <- corr_gauss_matrixC(x1, x2, th)
t2 <- corr_gauss_matrix_armaC(x1, x2, th)
identical(t1, t2)
microbenchmark::microbenchmark(corr_gauss_matrixC(x1, x2, th), corr_gauss_matrix_armaC(x1, x2, th))

x3 <- matrix(runif(1e3*6), ncol=)
t3 <- corr_gauss_matrix_symC(x3, th)
t4 <- corr_gauss_matrix_sym_armaC(x3, th)
identical(t3, t4)
microbenchmark::microbenchmark(corr_gauss_matrix_symC(x3, th), corr_gauss_matrix_sym_armaC(x3, th), times=50)
