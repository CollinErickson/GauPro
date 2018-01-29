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
