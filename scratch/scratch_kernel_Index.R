n <- 5
for (i in 0:(n-1-1)) {
  for (j in (i+1):(n-1)) {
    k = (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1
    print(c(i, j, k))
  }
}

kk <- IndexKernel$new(D=5)
kk
kk$kone(1,1, s2=1)
# kk$kone(1,2, s2=1)
kk$k(1,1)
kk$k(1,2)

kk$p <- (1:10)/100
kmat <- outer(1:5, 1:5, Vectorize(kk$k))
kmat <- matrix(NA,5,5)
for (i in 1:5) {
  for (j in 1:5) {
    kmat[i,j] <- kk$k(i,j)
    #print(c(i,j, kk$k(i,j)))
  }
}
kmat

kk$param_optim_lower()
kk$param_optim_upper()
kk$param_optim_start()
kk$param_optim_start0()


X <- matrix(c(1,1,1))


