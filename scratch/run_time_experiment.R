flist <- c(TestFunctions::banana,
           function(x) {TestFunctions::beambending(x)*1e6},
           TestFunctions::piston,
           TestFunctions::borehole)
e1 <- comparer::ffexp$new(
  n = c(20, 40, 80, 120, 160, 240, 320, 400, 500),
  FD = data.frame(
    d = c(2, 3, 7, 8),
    fi = 1:length(flist)
  ),
  k = c('Gaussian', 'Matern32'),
  eval_func = function(n, d, fi, k) {
    # browser()
    cat("\n", n, d, fi, k, "\n")
    f <- flist[[fi]]
    t1 <- proc.time()
    X <- matrix(runif(d*n), ncol=d)
    expect_no_warning(Z <- f(X))
    Z <- Z + rnorm(length(Z), 0, range(Z)*1e-4)
    gp <- GauPro_kernel_model$new(X, Z, kernel=k)
    t2 <- proc.time()
    (t2 - t1)[3]
  }
)
e1
e1$run_all(15)
e1$run_for_time(120, 1) #, run_order='shuffle')
e1$run_all(run_order = 'random')
e1
e1$plot()
lm(runtime ~ n + d + k, data=e1$outcleandf) %>% summary
lm(runtime^(1/3) ~ n + d + k, data=e1$outcleandf) %>% summary
lm(log(runtime) ~ n + d + k, data=e1$outcleandf) %>% summary
e1$outcleandf %>%
  ggplot(aes(n, log(runtime), color=d, shape=k)) +
  geom_point()
e1$outcleandf %>%
  ggplot(aes(n, runtime^(1/3), color=d, shape=k)) +
  geom_point()
e1$outcleandf$runtime %>% summary



# Can't get it to replace last of line
printfunc <- function() {
  cat("abcdef")
  Sys.sleep(1)
  cat("\r")
  # cat("ghi\n")
  cat("ghi\033[K\n")
}
printfunc()
