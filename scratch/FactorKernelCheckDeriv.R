nlev <- 3
kk <- LatentFactorKernel$new(D=1, nlevels=nlev, xindex=1, latentdim=2)
kk$p
kmat <- outer(1:nlev, 1:nlev, Vectorize(kk$k))
kmat
# Xmat <- matrix(sample(1:nlev, 19,T))
Xmat <- matrix(1:nlev, ncol=1)
kk$dC_dparams(X=Xmat, nug=0)
C1 <- kk$C_dC_dparams(X=Xmat, nug=0, params=c(kk$p, kk$s2))$C
C2 <- kk$k(Xmat, Xmat,params=c(kk$p, kk$s2))
max(abs(c(C1 - C2)))

kkpars <- c(kk$p, log(kk$s2,10))
kk$C_dC_dparams(X=Xmat, nug=0, params=kkpars)$C
eps <- 1e-6
kkparsplus <- kkpars
kkparsplus[1] <- kkparsplus[1] + eps
kkparsminus <- kkpars
kkparsminus[1] <- kkparsminus[1] - eps
cdcplus <- kk$C_dC_dparams(X=Xmat, nug=0, params=kkparsplus)
cdcminus <- kk$C_dC_dparams(X=Xmat, nug=0, params=kkparsminus)
der1 <- (cdcplus$C - cdcminus$C) / (2*eps)
der2 <- .5*(cdcplus$dC_dparams + cdcminus$dC_dparams)[1,,]
plot(c(der1), c(der2)); abline(a=0,b=1, col=2)
plot(c(der1), c(der1) - c(der2))


# Check product error
d <- 2
n <- 20
nlev <- 2
Xmat <- matrix(runif(d*n), ncol=2)
Xmat[,2] <- sample(1:nlev, n, T)
k1 <- IgnoreIndsKernel$new(Matern32$new(D=1, s2_est=F), ignoreinds = 2)
k2 <- LatentFactorKernel$new(D=2, nlevels = nlev, xindex = 2, latentdim = 1, s2_est = F)
k <- k1*k2
kpar1 <- k$k1$kernel$beta
kpar2 <- k$k2$p
kpar <- c(kpar1, kpar2)
k_C_dC <- k$C_dC_dparams(X=Xmat, nug=0, params = kpar)
eps <- 1e-6
for (i in 1:length(kpar)){
  kkparsplus <- kpar
  kkparsplus[i] <- kkparsplus[i] + eps
  kkparsminus <- kpar
  kkparsminus[i] <- kkparsminus[i] - eps
  debugonce(k$k2$dC_dparams)
  cdcplus <- k$C_dC_dparams(X=Xmat, nug=0, params=kkparsplus)
  cdcminus <- k$C_dC_dparams(X=Xmat, nug=0, params=kkparsminus)
  der1 <- (cdcplus$C - cdcminus$C) / (2*eps)
  der2 <- .5*(cdcplus$dC_dparams + cdcminus$dC_dparams)[i,,]
  print(i); print(summary(c(der1 - der2)))
  plot(c(der1), c(der2)); abline(a=0,b=1, col=2)
  plot(c(der1), c(der1) - c(der2))
}
