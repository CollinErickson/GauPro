set.seed(0)
n <- 20
x <- lhs::maximinLHS(n=n, k=2)
y <- TestFunctions::banana(x)
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian, nug=1e-8, nug.est = F)
# xx <- matrix(runif(2e5),ncol=2)
# gp$pred(XX = xx, se=T)

gp$plot2D()
x1 <- c(.6,.83)
x2 <- c(.58,.83)
set.seed(1)
print(gp$pred(x1,se=T)$s2) # 0.1749093
print(gp$pred(x2,se=T)$s2) # 0.1708879
print(gp$pred_var_after_adding_points(add_points = x1, pred_points = x2)) # 0.02812571
print(gp$pred_var_reduction(add_point = x1, pred_points = x2))

gp$update(Xnew = x1, Znew = TestFunctions::banana(x1), no_update = TRUE)
print(gp$pred(x1,se=T)$s2) # 0.1749093
print(gp$pred(x2,se=T)$s2) # 0.1708879

# Check with multiple points
x3 <- matrix(c(.61,.83,.62,.83), 2,2, byrow = T)
x4 <- matrix(c(.61,.84,.62,.84), 2,2, byrow = T)
print(gp$pred(x4, se=T)$s2)
print(gp$pred_var_after_adding_points(add_points = x3, pred_points = x4))
# Actually add them
gp$update(Xnew = x3, Znew = TestFunctions::banana(x3), no_update = TRUE)
print(gp$pred(x4,se=T)$s2)


# Create two more
set.seed(2)
x5 <- matrix(runif(2*3), ncol=2)
x6 <- matrix(runif(1e3*2), ncol=2)
print(gp$pred(x5, se=T)$s2)
print(gp$pred_var_after_adding_points(add_points = x5, pred_points = x6))
# Time
print(microbenchmark::microbenchmark({gp$pred_var_after_adding_points(add_points = x5, pred_points = x6)}, times=100))
# profvis
# pv1 <- profvis::profvis(gp$pred_var_after_adding_points(add_points = x5, pred_points = x6))
# pv1 <- profvis::profvis(replicate(100,gp$pred_var_after_adding_points(add_points = x5, pred_points = x6)))

# Time
print(microbenchmark::microbenchmark({gp$pred_var_reduction(add_point = x5[1,], pred_points = x6)}, times=100))
# pv1 <- profvis::profvis(replicate(100,gp$pred_var_reduction(add_point = x5[1,], pred_points = x6)))

print(microbenchmark::microbenchmark({
  gp$update(Xnew = x5, Znew = TestFunctions::banana(x5), no_update = TRUE)
  (gp$pred(x6,se=T)$s2)
}, times=1))
print((gp$pred(x6,se=T)$s2))







set.seed(0)
n <- 30
x <- lhs::maximinLHS(n=n, k=2)
y <- TestFunctions::banana(x)
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian, nug=1e-8, nug.est = F)
xx <- matrix(runif(2e5),ncol=2)
x1 <- matrix(c(.3,.4), ncol=2)
gp$pred_var_after_adding_points(add_points = x1, pred_points = xx)
pv5 <- profvis::profvis(gp$pred_var_after_adding_points(add_points = x1, pred_points = xx))

# Had made 3 different methods, deleted 2 slowest, kept fastest
#method1 <- 1
system.time(t1 <- gp$pred_var_after_adding_points(add_points = x1, pred_points = xx))
#method1 <- 2
#system.time(t2 <- gp$pred_var_after_adding_points(add_points = x1, pred_points = xx))
#method1 <- 3
#system.time(t3 <- gp$pred_var_after_adding_points(add_points = x1, pred_points = xx))

# Check matrix vec
x2 <- c(.33,.43)
gp$pred_var_after_adding_points(add_points = x2, pred_points = xx)
gp$pred_var_after_adding_points(add_points = x2, pred_points = c(.5,.6))

gp$pred_var_reduction(add_point = x2, pred_points = x1)
gp$pred_var_reduction(add_point = x2, pred_points = c(.5,.6))
gp$pred_var_reduction(add_point = x2, pred_points = matrix(c(.5,.6), ncol=2))
gp$pred_var_reduction(add_point = x2, pred_points = xx)
apply(xx, 1, function(xi)gp$pred_var_reduction(add_point = x2, pred_points = xi))


# Test pred_var_reductions, vectorization of pred_var_reduction
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian, nug=1e-8, nug.est = F)
xx <- matrix(runif(2e3),ncol=2)
yy <- matrix(runif(200), ncol=2)
t1 <- gp$pred_var_reductions(add_points = yy, pred_points = xx)
t2 <- apply(yy, 1, function(yi) gp$pred_var_reduction(add_point = yi, pred_points = xx))
str(t1)
str(t2)
summary(c(t1 - t2))
hist(c(t1-t2), breaks=40)
# Check (7,1) of t1-t2 since it is .068 diff
gp$pred_var_reduction(add_point = yy[1,,drop=T], pred_points = xx[7,,drop=F])
gp$pred_var_reductions(add_points = yy[1,,drop=F], pred_points = xx[7,,drop=F])
i <- sample(1:100,1);j <- sample(1:100,1); gp$pred_var_reduction(add_point = yy[j,,drop=T], pred_points = xx[i,,drop=F]) - gp$pred_var_reductions(add_points = yy[j,,drop=F], pred_points = xx[i,,drop=F])

# Use smaller samples
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian, nug=1e-8, nug.est = F)
xx <- matrix(runif(2*3),ncol=2)
yy <- matrix(runif(2*2), ncol=2)
t1 <- gp$pred_var_reductions(add_points = yy, pred_points = xx)
t2 <- apply(yy, 1, function(yi) gp$pred_var_reduction(add_point = yi, pred_points = xx))
str(t1)
str(t2)
summary(c(t1 - t2))
# Try to debugo
debugonce(gp$pred_var_reduction); apply(yy, 1, function(yi) gp$pred_var_reduction(add_point = yi, pred_points = xx))
debugonce(gp$pred_var_reductions); gp$pred_var_reductions(add_points = yy, pred_points = xx)

# Fixed it!!!
# Now benchmark
microbenchmark::microbenchmark(
  t1 <- gp$pred_var_reductions(add_points = yy, pred_points = xx),
  t2 <- apply(yy, 1, function(yi) gp$pred_var_reduction(add_point = yi, pred_points = xx)))
# This is 10x faster
# Unit: milliseconds
# expr       min        lq      mean
# t1 <- gp$pred_var_reductions(add_points = yy, pred_points = xx)  4.664792  5.092234  11.08233
# t2 <- apply(yy, 1, function(yi) gp$pred_var_reduction(add_point = yi,      pred_points = xx)) 85.510294 87.586110 104.24953
# median         uq      max neval cld
# 5.981516   7.112744 170.2368   100  a
# 89.599399 103.002694 290.8661   100   b



# Check pred_var_after_adding_points_sep
set.seed(0)
n <- 29
x <- lhs::maximinLHS(n=n, k=2)
y <- TestFunctions::banana(x)
xx <- lhs::randomLHS(n=30, k = 2)
yy <- lhs::randomLHS(n=1e4, k = 2)
gp <- GauPro_kernel_model$new(X=x, Z=y, kernel=Gaussian, nug=1e-8, nug.est = F)
t1 <- apply(xx, 1, function(xi) gp$pred_var_after_adding_points(add_points = xi, pred_points = yy))
t2 <- gp$pred_var_after_adding_points_sep(add_points = xx, pred_points = yy)
summary(c(t1-t2))
microbenchmark::microbenchmark(
  apply(xx, 1, function(xi) gp$pred_var_after_adding_points(add_points = xi, pred_points = yy))
  , gp$pred_var_after_adding_points_sep(add_points = xx, pred_points = yy)
)
# 3.5x faster use sep version

# Checking growth, its linear in add_points and pred_points
# ms <- c(100, 300,1e3, 2e3, 3e3, 5e3, 6e3, 7e3, 8e3, 1e4, 2e4)
ms <- c(20,40,50,60,70,80,90,1e2)
xx2 <- matrix(runif(2000), ncol=1)
ts <- sapply(
  ms,
  function(m) {
    xm <- lhs::randomLHS(n = m, k = 2)
    gp$update(Xall = xm, Zall = TestFunctions::banana(xm))
    system.time({
      # Check growth in add_points
      # gp$pred_var_after_adding_points_sep(add_points = xm, pred_points = xx)
      # Check growth in pred_points
      # gp$pred_var_after_adding_points_sep(add_points = xx, pred_points = xm)
      # Check growth in design points
      gp$pred_var_after_adding_points_sep(add_points = xx, pred_points = xm)
    })[3]
  }
)
plot(ms, ts)
