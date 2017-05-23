set.seed(1)
n <- 20
d <- 1
f <- function(x) sin(2*pi*x)^.8
f <- function(x) abs(sin(2*pi*x))^.8 * sign(sin(2*pi*x)) * x^4 + rnorm(n, 0, .001)
X <- lhs::maximinLHS(n, d) #matrix(runif(n*d), ncol=2)
Z <- f(X)
gp <- GauPro(X, Z)#, nug=1e-8, nug.est=F)

E <- gp$Kinv[1:(n-1), 1:(n-1)]
B <- gp$K[n, 1:(n-1)]
G <- gp$Kinv[n, 1:(n-1)]
Ainv1 <- E + E %*% outer(B, G) / (1 - sum(B*G))

Ainv2 <- solve(gp$K[-n, -n])

microbenchmark::microbenchmark(solve(gp$K[-n, -n]),
                               {
                                 E <- gp$Kinv[1:(n-1), 1:(n-1)]
                                 B <- gp$K[n, 1:(n-1)]
                                 G <- gp$Kinv[n, 1:(n-1)]
                                 E + E %*% outer(B, G) / (1 - sum(B*G))
                               },
                               {
                                 gp$Kinv[-n, -n] %*% (diag(n-1) + outer(gp$K[n, 1:(n-1)], gp$Kinv[n, 1:(n-1)]) / (1 - sum(gp$K[n, 1:(n-1)]*gp$Kinv[n, 1:(n-1)])))
                               },
                               {
                                 gp$Kinv[-n, -n]+ gp$Kinv[-n, -n] %*% gp$K[n, 1:(n-1)] %*% gp$Kinv[n, 1:(n-1)] / (1 - sum(gp$K[n, 1:(n-1)]*gp$Kinv[n, 1:(n-1)]))
                               },
                               {
                                 gp$Kinv[-n, -n]+ gp$Kinv[-n, -n] %*% gp$K[n, -n] %*% gp$Kinv[n, -n] / (1 - sum(gp$K[n, -n]*gp$Kinv[n, -n]))
                               },
                               times = 100)

yhati <- function(i, gp) {#browser()
  #n <- nrow(gp$X)
  Kinvi <- gp$Kinv[-i, -i]+ gp$Kinv[-i, -i] %*% gp$K[i, -i] %*% gp$Kinv[i, -i] / (1 - sum(gp$K[i, -i]*gp$Kinv[i, -i]))
  mn <- gp$mu_hat + sum(gp$K[i, -i] %*% Kinvi * (gp$Z[-i] - gp$mu_hat))
  vr <- pmax(1e-16, gp$K[i, i] - gp$K[i, -i] %*% Kinvi %*% gp$K[-i, i])
  c(mn, vr)
}
yhats <- sapply(1:n, function(i){yhati(i, gp)})
#plot(yhats, Z)
ContourFunctions::cf(X, Z-yhats[1,])
zhats <- (yhats[1,] - Z) / sqrt(yhats[2,])
abszhats <- abs(zhats)
gpz <- GauPro(X, abszhats, nug.est = F, nug = 1e-8, theta = gp$theta, param.est=F)

predse <- function(XX) {#browser()
  pr <- gp$pred(XX, se.fit = T)
  pr_z <- pmax(1e-4, gpz$pred(XX))
  pr_se <- pr$se * pr_z
  #pr$s2 <- pr$s2 * pr_z^2
  pr_se
}
ContourFunctions::cf(gp$pred)
ContourFunctions::cf(function(x)gp$predict(x, se=T)$se, batchmax=Inf, pts=X)
ContourFunctions::cf(function(x)predse(x), batchmax=Inf, pts=X)
XX <- lhs::randomLHS(1e4, d)
ZZ <- f(XX)
ZZpred <- gp$pred(XX)
ZZse1 <- gp$pred(XX, se=T)$se
ZZse2 <- predse(XX)
plot(abs(ZZpred - ZZ), ZZse1);abline(a=0,b=1,col=2)
plot(abs(ZZpred - ZZ), ZZse2);abline(a=0,b=1,col=2)
hist((ZZpred - ZZ) / ZZse1)
hist((ZZpred - ZZ) / ZZse2)

curve(gp$pred(x))
points(X, Z)
curve(gp$pred(x) + gp$pred(x,se=T)$se, add=T, col=2)
curve(gp$pred(x) - gp$pred(x,se=T)$se, add=T, col=2)
curve(gp$pred(x) + predse(x), add=T, col=3)
curve(gp$pred(x) - predse(x), add=T, col=3)
