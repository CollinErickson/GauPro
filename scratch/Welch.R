tf <- function(x) {
  5*x[12]/(1+x[1]) + 5*(x[4]-x[20])^2 + x[5] + 40*x[19]^3 - 5*x[19] +
    .05*x[2] + .08*x[3] -.03*x[6] + .03*x[7] -.09*x[9] -.01*x[10] -.07*x[11] +
    .25*x[13]^2 -.04*x[14] +.06*x[15] -.01*x[17] -.03*x[18]
}
X <- lhs::maximinLHS(n=100,k=20) -.5
Z <- apply(X, 1, tf)
gp <- GauPr_Gauss$new(X=X, Z=Z, theta_map=rep(1,20), nug.est=F, nug=1e-8)


dev0 <- gp$deviance()
devs <- rep(Inf, 20)
llh0 <- gp$loglikelihood()
llhs <- rep(-Inf, 20)
avail.ind <- 1:20
for (i in 1:5) {
  for (d in avail.ind) {
    gpd <- gp$clone()
    if (sum(gpd$theta_map[d] == gpd$theta_map) == 1) { # already unique
      # pass
    } else { # not unique
      #gpd$nug <- 1
      gpd$theta_short <- c(gpd$theta_short, gpd$theta_short[gpd$theta_map[d]])
      gpd$theta_length <- gpd$theta_length + 1
      gpd$theta_map[d] <- max(gpd$theta_map) + 1
      gpd$update(restarts=20)
      devs[d] <- gpd$deviance()
      llhs[d] <- gpd$loglikelihood()
    }
    rm(gpd)
  }
  ind <- which.max(llhs)
  #ind <- which.min(devs)
  llh <- llhs[ind]
  dev <- devs[ind]
  inc2 <- 2*(llh-llh0)
  print(paste('best is',ind, 'llh is', llh, '2inc is', inc2))
  gp$theta_short <- c(gp$theta_short, gp$theta_short[gp$theta_map[ind]])
  gp$theta_length <- gp$theta_length + 1
  gp$theta_map[ind] <- max(gp$theta_map) + 1
  gp$update()
  dev0 <- dev
  llh0 <- llh
  devs <- rep(Inf, 20)
  llhs <- rep(-Inf, 20)
  avail.ind <- setdiff(avail.ind, ind)
}
