# Scratch to test getting bernoulli working
n <- 50
d <- 2
logistic <- function(x) 1/(1+exp(-x))
dlogistic <- function(x) logistic(x) * (1 - logistic(x))
logit <- function(p) log(p/(1-p))
# Unit vector
ei <- function(i) rep(1,n) * (i == (1:n))

X <- matrix(runif(n*d), nrow=n)
# w_true <- apply(X, 1, function(x) (2*(x[1]-.5)) + 2*sin(x[2]*pi) - 1)
w_true <- apply(X, 1, function(x) (4*(x[1]-.5)))
logistic_w_true <- logistic(w_true)
y <- as.integer(runif(n) < logistic_w_true)
cbind(X, w_true, logistic_w_true, y)
cor(cbind(X, w_true, logistic_w_true, y))
pairs(cbind(X, w_true, logistic_w_true, y))

w <- rnorm(n, sd=1.7)
w
kern <- GauPro::k_Gaussian(s2=1.7^2, s2_est=F, D=d)
kern$plot()
gp <- gpkm(kern=kern, X, w)
gp$sample(X)
gp$plotmarginal()
w_llh <- function(u) {
  -0.5*gp$D*log(2*pi) - 0.5*determinant(gp$K, logarithm = T)$modulus[1] +
    -0.5*sum((u - gp$mu_hatX) * solve(gp$K, u - gp$mu_hatX))
}
w_llh(w)
w_llh_gr <- function(u) {
  -0.5 * (2 * solve(gp$K, (u - gp$mu_hatX)))[,1]
}
# Check grad, is right
list((w_llh(w + ei(50)*eps) - w_llh(w))/eps, w_llh_gr(w))


y_llh <- function(u) {
  sum(log(ifelse(y > 0.5, logistic(u), 1-logistic(u))))
}
y_llh(w)
y_llh(w_true)
y_llh_gr <- function(u) {
  (ifelse(y > 0.5,
              dlogistic(u) / logistic(u),
             -dlogistic(u)/ (1 - logistic(u))))
}
# Check grad, it's right
list((y_llh(w + ei(3)*eps) - y_llh(w))/eps, y_llh_gr(w))
# Check that y_llh optim works with grad, it does, can't use default method
optim(w, y_llh, y_llh_gr, control=list(fnscale=-1), method='BFGS')

predict_y <- function(x) {
  w_pred <- gp$predict(x)
  logistic(w_pred)
}
plot(predict_y(X), w_true)
# Loop over two optimization steps
for (i in 1:30) {
  cat('step', i, mean(abs(w - w_true)), "\t", mean(abs(y - logistic(w))), '\t', mean(w), "\n")
  # Optimization 1: Optimize w, keep GP params constant
  optim_w_func <- function(u) {
    w_llh(u) + y_llh(u)
  }
  optim_w_func_gr <- function(u) {
    w_llh_gr(u) + y_llh_gr(u)
  }
  cat('\tOpt 1', optim_w_func(w), '\t', w_llh(w), '\t', y_llh(w), "\n")
  # w <- optim(w, optim_w_func, control=list(fnscale=-1))$par
  w <- optim(w, optim_w_func, optim_w_func_gr, method='BFGS', control=list(fnscale=-1))$par
  cat('\tOpt 1', optim_w_func(w), '\t', w_llh(w), '\t', y_llh(w), "\n")

  # Optimization 2: Optimize GP params, keep w constant
  cat('\t\tOpt 2', gp$loglikelihood(), "\n")
  gp$update(Zall = w)
  cat('\t\tOpt 2', gp$loglikelihood(), "\n")
}
plot(w, w_true)
plot(X[,1], y)
points(X[,1], logistic(w_true), col=2)
# This looks wrong
points(X[,1], logistic(w), col=3)
points(X[,1], logistic(gp$pred(X)), col=4)
