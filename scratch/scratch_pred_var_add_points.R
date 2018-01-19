set.seed(0)
n <- 50
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
x5 <- matrix(c(.61,.83,.62,.83), 2,2, byrow = T) - .2
x6 <- matrix(c(.61,.84,.62,.84), 2,2, byrow = T) - .2
print(gp$pred(x5, se=T)$s2)
print(gp$pred_var_after_adding_points(add_points = x5, pred_points = x6))
# Time
print(microbenchmark::microbenchmark({gp$pred_var_after_adding_points(add_points = x5, pred_points = x6)}, times=1))

print(microbenchmark::microbenchmark({
  gp$update(Xnew = x5, Znew = TestFunctions::banana(x5), no_update = TRUE)
  (gp$pred(x6,se=T)$s2)
}, times=1)
