
kk <- LatentFactorKernel$new(D=1, nlevels=5, xindex=1, latentdim=2)
kk$p
kmat <- outer(1:5, 1:5, Vectorize(kk$k))
kmat
kk$dC_dparams(X=matrix(1:5, ncol=1), nug=0)
kk$C_dC_dparams(X=matrix(1:5, ncol=1), nug=0, params=c(kk$p, kk$s2))$C


kkpars <- c(kk$p, log(kk$s2,10))
kk$C_dC_dparams(X=matrix(1:5, ncol=1), nug=0, params=kkpars)$C
eps <- 1e-1
kkparsplus <- kkpars
kkparsplus[1] <- kkparsplus[1] + eps
kkparsminus <- kkpars
kkparsminus[1] <- kkparsminus[1] - eps
cdcplus <- kk$C_dC_dparams(X=matrix(1:5, ncol=1), nug=0, params=kkparsplus)
cdcminus <- kk$C_dC_dparams(X=matrix(1:5, ncol=1), nug=0, params=kkparsminus)
(cdcplus$C - cdcminus$C) / (2*eps)
.5*(cdcplus$dC_dparams + cdcminus$dC_dparams)[1,,]
