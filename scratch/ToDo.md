# GauPro to do

* Factor kernel: example

* Ordered factor kernel: example

* Ignore input kernel: example

* Calculate EI: nopt, test, doc

* EI should use mu instead of Z for noise. Should it be sd of mean, not pred?
See http://krasserm.github.io/2018/03/21/bayesian-optimization/.

* max_qEI: test, doc, exact solution. Change CL to CLpred

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

* 4 inputs, 2 latentfac, .1 noise. predictions have near zero s2

* factor trend, LM ignore inds

* diagnostics

* transform X, Z

* penalties in optim, on trend/kernels

* Clean up trend. b_est/m_est, jitter, aren't used

* Check which indices are non-factor. If any overlap, give warning.

* All restarts had error: doesn't say what error is, so sometimes it's
an error the user should know. Like when I test useCM.

* dEI_dx

* maxEI prob doesn't work because the differences are tiny. Log scale?

* maxqEI doesn't spread out much

* Knowledge gradient: multiple starts for optim

* Change plot for factor inputs

* When giving in formula/data
  * fix EI
  * fix plot for 1D, 2D, marginal, marginal random
  * fix
  * convert ordered factor
  * add to doc
  * convert back for formula, need for EI
  * make sure exp(), log(), *, : work for fit and pred. Should * be allowed?
  * message what chosen kernel is

* maxqEI: if factor in kernel and no mopar, give warning.

* FactorKernel and OrderedFK: change correlation, maybe just 1.95 instead of 2

* mixopt for qEI: picks same point multiple times b/c of mean uncertainty

* Should EI be zero for already evaluated points? Mean can have high uncertainty
when nug>0, so it can pick already eval points.

* Resubmit to CRAN by 11/28. After mixopt update.

* Plot doesn't work for gp4$kernel

* gpkm summary: diagnostics from LOO, check LM for ideas 
