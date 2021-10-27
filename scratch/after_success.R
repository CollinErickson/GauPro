timestamp()
cat("RUNNING after_success.R ...\n")
timestart <- Sys.time()

# borehole <- function(x) {
#   rw <- x[, 1] * (0.15 - 0.05) + 0.05
#   r <-  x[, 2] * (50000 - 100) + 100
#   Tu <- x[, 3] * (115600 - 63070) + 63070
#   Hu <- x[, 4] * (1110 - 990) + 990
#   Tl <- x[, 5] * (116 - 63.1) + 63.1
#   Hl <- x[, 6] * (820 - 700) + 700
#   L <-  x[, 7] * (1680 - 1120) + 1120
#   Kw <- x[, 8] * (12045 - 9855) + 9855
#
#   m1 <- 2 * pi * Tu * (Hu - Hl)
#   m2 <- log(r / rw)
#   m3 <- 1 + 2 * L * Tu / (m2 * rw ^ 2 * Kw) + Tu / Tl
#   return(m1 / m2 / m3)
# }


d = 7
testf <- TestFunctions::piston #borehole #borehole

N <- 50
Npred <- 1000


X <- matrix(runif(N*d), N, d)
Y = testf(X)

Xp <- matrix(runif(Npred*d), Npred, d)
Yp = testf(Xp)

require("GauPro")


fittimestart <- Sys.time()
SG = GauPro::GauPro_kernel_model$new(X=X, Z=Y, kernel=GauPro::Gaussian) #do one final parameter estimation
fittimeend <- Sys.time()
# SG = CGGPcreate(d=d, batchsize=201)
# Y = testf(SG$design)
# SG = CGGPfit(SG, Y)
# cat("Now doing Bayesian\n")

# # Changing this so it doesn't call fit 40 times with no param updates
# for (ic in (1:round((N-201)/200))[1:9]) {
#   cat(ic, " ")
#   SG=CGGPappend(SG,200, "MAP") #add 200 points to the design based on thetahat
#   Y = testf(SG$design)
#   if( ic< 10){  #eventually we stop estimating theta because it takes awhile and the estimates dont change that much
#     SG = CGGPfit(SG,Y) #estimate the parameter (SG structure is important)
#   }
# }
# SG=CGGPappend(SG, 200*length((1:round((N-201)/200))[-(1:9)]), "MAP")

# cat("\n")
# Y = testf(SG$design)
# timelastlogthetaMLEstart <- Sys.time()
# SG = CGGPfit(SG, Y) #do one final parameter estimation
# timelastlogthetaMLEend <- Sys.time()

timepredstart <- Sys.time()
GP = predict(SG, Xp, se.fit=T) #build a full emulator
timepredend <- Sys.time()

RMSE <- sqrt(mean(((Yp-GP$mean)^2)))  #prediction should be much better
meanscore <- mean((Yp-GP$mean)^2/GP$s2+log(GP$s2)) #score should be much better
meancoverage <- mean((Yp<= GP$mean+1.96*(GP$se))&(Yp>= GP$mean-1.96*(GP$se)))  #coverage should be closer to 95 %

cat("RMSE is         ", RMSE, "\n")
cat("R-sq is         ", 1-RMSE^2/mean((mean(Y)-Yp)^2), "\n")
cat("Mean score is   ", meanscore, "\n")
cat("coverage is     ", meancoverage, "\n")

# Don't count plotting in run time
timeend <- Sys.time()
cat("Total run   time is:", capture.output(timeend - timestart), '\n')
cat("Prediction  time is:", capture.output(timepredend - timepredstart), '\n')
cat("Fit time is:", capture.output(fittimeend - fittimestart), '\n')

# if (T) { # Can Travis just skip this?
#   di <- sample(1:nrow(SG$design), 100)
#   Y0pred <- CGGPpred(SG$design[di,],CGGP=SG) #,Y,logtheta=logthetaest)
#   plot(Yp, GP$mean, ylim=c(min(GP$mean, Y0pred$m),max(GP$mean, Y0pred$m))); points(Y[di], Y0pred$m,col=3,pch=2); abline(a=0,b=1,col=2)
#   # Now plot with bars
#   #plot(Yp, GP$mean , ylim=c(min(GP$mean, Y0pred$m),max(GP$mean, Y0pred$m)),pch=19)#; points(Y[di], Y0pred$m,col=3,pch=2); abline(a=0,b=1,col=2)
#   plot(Yp, Yp-GP$mean , ylim=max(sqrt(GP$var))*c(-2,2))#c(min(-2GP$mean, Y0pred$m),max(GP$mean, Y0pred$m)),pch=19,col='white')#; points(Y[di], Y0pred$m,col=3,pch=2); abline(a=0,b=1,col=2)
#   points(Yp, 0*GP$mean + 2*sqrt(GP$var), col=4, pch=19)#; points(Y[di], Y0pred$m,col=3,pch=2); abline(a=0,b=1,col=2)
#   points(Yp, 0*GP$mean - 2*sqrt(GP$var), col=5)#; points(Y[di], Y0pred$m,col=3,pch=2); abline(a=0,b=1,col=2)
#   errmax <- max(sqrt(GP$var), abs(GP$mean - Yp))
#   plot(GP$mean-Yp, sqrt(GP$var), xlim=errmax*c(-1,1), ylim=c(0,errmax))#;abline(a=0,b=1,col=2)
#   polygon(1.1*errmax*c(0,-2,2),1.1*errmax*c(0,1,1), col=3, density=10, angle=135)
#   polygon(1.1*errmax*c(0,-1,1),1.1*errmax*c(0,1,1), col=2, density=30)
#   points(GP$mean-Yp, sqrt(GP$var), xlim=errmax*c(-1,1), ylim=c(0,errmax))
# }

# cat(capture.output(Sys.time() - timestart), '\n')
cat("... FINISHED after_success.R\n")
timestamp()

