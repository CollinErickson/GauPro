# GauPro to do

* Factor kernel: example

* Ordered factor kernel: example

* Ignore input kernel: example

* Take in data/formula input in addition to matrix

* Calculate EI: nopt, test, doc

* EI should use mu instead of Z for noise. Should it be sd of mean, not pred?
See http://krasserm.github.io/2018/03/21/bayesian-optimization/.

* max_qEI: test, doc, exact solution. Change CL to CLpred

* With many points, just estimate param from subset of data, but use all data
at end for K/Kinv so predictions are good.

* Change print, summary

* Improve kernel model plot

* Give kernel in as string. If missing, give better error or just use Mat52.

* Reduce L-BFGS-B tolerance? Need to change tests too.

* Add documentation for kernels, esp. factor ones

* Better handling of ignore input when using factors?

* Transform inputs

* Change S3 plot to R6

* plot2D option to do se instead of mean, or side by side

* use self$kernel$k on x instead of assuming it is 1

* Fix crashing on laptop. Changed gradfuncarray from Rcpp to R, but something
in optimRestart is now breaking instead.

* optim NaN starting value (found on 1D doing EI)

* Speed up m52 grad, Cubic k/grad

* 3 kernel product. prod/sum need to get s2_est from children

* no param est gave error, dC_dparam can't be calculated

* Add test with repeated X. Add test with big nugget and make sure deriv still matches.

* estimate_run_time()

* test factor kernels, check grad

* data frame input: everywhere, predict, etc. Convert to matrix.

* 4 inputs, 2 latentfac, .1 noise. predictions have near zero s2

* speed up latentfactor (dC in Rcpp)

* factor trend

* diagnostics

* transform X, Z

* penalties in optim, on trend/kernels

* Clean up trend. b_est/m_est, jitter, aren't used

* kernel print

* track optim: track deviance

* Check which indices are non-factor. If any overlap, give warning.

* why does predict return a matrix? Just make it a vector.
