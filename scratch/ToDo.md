# GauPro to do

* Factor kernel: example

* Ordered factor kernel: example

* Ignore input kernel: example

* Calculate EI: nopt, test, doc

* EI should use mu instead of Z for noise. Should it be sd of mean, not pred?
See http://krasserm.github.io/2018/03/21/bayesian-optimization/.

* With many points, just estimate param from subset of data, but use all data
at end for K/Kinv so predictions are good.

* Improve kernel model plot

* Reduce L-BFGS-B tolerance? Need to change tests too.

* Add documentation for kernels, esp. factor ones

* Better handling of ignore input when using factors?

* Transform inputs

* plot2D option to do se instead of mean, or side by side

* use self$kernel$k on x instead of assuming it is 1

* optim NaN starting value (found on 1D doing EI)

* Speed up m52 grad, Cubic k/grad

* no param est gave error, dC_dparam can't be calculated

* Add test with repeated X. Add test with big nugget and make sure deriv still matches.

* estimate_run_time()

* test factor kernels, check grad

* factor trend, LM ignore inds

* transform X, Z

* penalties in optim, on trend/kernels

* Clean up trend. b_est/m_est, jitter, aren't used

* Check which indices are non-factor. If any overlap, give warning.

* All restarts had error: doesn't say what error is, so sometimes it's
an error the user should know. Like when I test useCM.

* dEI_dx

* maxEI prob doesn't work because the differences are tiny. Log scale?

* Knowledge gradient: multiple starts for optim

* When giving in formula/data
  * fix plot for 1D, 2D, marginal, marginal random
  * convert ordered factor
  * add to doc
  * message what chosen kernel is

* OrderedFK and Latent: change correlation, maybe just 1.95 instead of 2

* qEI with mixopt:
    * picks same point multiple times b/c of mean uncertainty
    * error converting mopar back and forth. mopar_converted?
    * Slow to convert to df each time. Convert mopar to mopar_converted.
    * doesn't spread out much. Better with AugEI or CorEI?
    * test, doc.
    * CL or pred?

* Plot doesn't work for gp4$kernel

* Plot between marginal and marginalrandom. Vary along 10 points, partial range.
Hard to get factors right.

* Corrected EI

* Augmented EI: add, grad, test

* Is EI minimize arg working properly?

* Optimize any function. Avoid reimplementing maxEI, maxAugEI, max___.
  * Can use gr or fngr
  * Can do matrix eval
  * mopar or mopar converted

* Large Z variance is bad. Extend range of s2. Warning to rescale. Tell kernel
that initial optim values should be large.

* If after par optim, par are at lower/upper, give warning
