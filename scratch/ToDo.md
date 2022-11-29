# GauPro to do

* EI: nopt, test, doc
  * Use t-dist
  * dEI_dx

* With many points, just estimate param from subset of data, but use all data
at end for K/Kinv so predictions are good.

* Reduce L-BFGS-B tolerance? Need to change tests too.

* Add documentation for kernels, esp. factor ones

* Better handling of ignore input when using factors?

* Transform inputs

* optim NaN starting value (found on 1D doing EI)

* Speed up m52 grad, Cubic k/grad

* no param est gave error, dC_dparam can't be calculated

* Add test with repeated X. Add test with big nugget and make sure deriv still matches.

* estimate_run_time()

* factor trend, LM ignore inds

* transform X, Z

* penalties in optim, on trend/kernels

* Clean up trend. b_est/m_est, jitter, aren't used

* Check which indices are non-factor. If any overlap, give warning.

* All restarts had error: doesn't say what error is, so sometimes it's
an error the user should know. Like when I test useCM.

* Knowledge gradient: multiple starts for optim

* When giving in formula/data
  * fix plot for 1D, 2D, marginal, marginal random
  * convert ordered factor
  * add to doc
  * message what chosen kernel is

* qEI with mixopt:
    * picks same point multiple times b/c of mean uncertainty
    * error converting mopar back and forth. mopar_converted?
    * Slow to convert to df each time. Convert mopar to mopar_converted.
    * doesn't spread out much. Better with AugEI or CorEI?
    * test, doc.
    * CL or pred?

* Corrected EI

* Augmented EI: grad, test, check minimize

* Optimize any function. Avoid reimplementing maxEI, maxAugEI, max___.
  * Can use gr or fngr
  * Can do matrix eval
  * mopar or mopar converted

* Large Z variance is bad. Extend range of s2. Warning to rescale. Tell kernel
that initial optim values should be large.

* If after par optim, par are at lower/upper, give warning

* Add deprecated for old model: lifecycle::

* Add test with factor kernel in multiple kernels, or else make it give error
* Make sure no overlap with factor and cts indexes

* Use blockcd, not coordesc

* Submit bug for lifecycle depr inside R6
