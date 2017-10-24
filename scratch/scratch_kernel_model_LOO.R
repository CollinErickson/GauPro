
# Check banana again
n <- 40
d <- 2
f1 <- TestFunctions::banana#function(x) {abs(sin(2*pi*x[1]))}
X1 <- MaxPro::MaxProLHD(n=n,p=d, itermax = 10)$Design#t(lhs::maximinLHS(n=d,k=n)) #
# X1 <- matrix(runif(n*d),n,d)
Z1 <- apply(X1,1,f1) * 9.3 # + rnorm(n, 0, 1e-3)
gp <- GauPro_kernel_model_LOO$new(X=X1, Z=Z1, kernel=Exponential)
ContourFunctions::cf(gp$predict, pts=X1)
gp_noLOO <- gp$clone(deep = T); gp_noLOO$use_LOO <- FALSE
ContourFunctions::cf(gp_noLOO$predict, pts=X1)
ContourFunctions::cf(gp$tmod$predict, pts=X1)
# Plot se's
ContourFunctions::cf(function(x) gp$pred(x, se=T)$se, pts=X1, batchmax=Inf)
ContourFunctions::cf(function(x) gp_noLOO$pred(x, se=T)$se, pts=X1, batchmax=Inf)
# See predicted t values
nn <- 1e3
XX <- matrix(runif(nn*d),nn,d)
ZZ <- apply(XX, 1, f1) * 9.3
ploo <- gp$pred(XX, se=T)
tloo <- (ploo$mean - ZZ) / ploo$se
summary(tloo)
pnoloo <- gp_noLOO$pred(XX, se=T)
tnoloo <- (pnoloo$mean - ZZ) / pnoloo$se
summary(tnoloo)
stripchart(list(LOO=tloo,no_LOO=tnoloo))
ContourFunctions::cf(function(x){z <- apply(x, 1, f1)*9.3; ploo <- gp$pred(x, se=T);tloo <- (ploo$mean - z) / ploo$se;abs(tloo)}, pts=X1, batchmax=Inf)
ContourFunctions::cf(function(x){z <- apply(x, 1, f1)*9.3; ploo <- gp_noLOO$pred(x, se=T);tloo <- (ploo$mean - z) / ploo$se;abs(tloo)}, pts=X1, batchmax=Inf)
qqnorm(tloo)
qqline(tloo)
qqnorm(tnoloo)
qqline(tnoloo)
plot(tloo, gp$tmod$predict(XX))
