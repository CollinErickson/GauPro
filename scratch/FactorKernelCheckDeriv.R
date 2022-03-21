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
