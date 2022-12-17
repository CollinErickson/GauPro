# GauPro to do

* EI: nopt, test, doc
  * Use t-dist
  * Use dEI_dx in maxEI

* With many points, just estimate param from subset of data, but use all data
at end for K/Kinv so predictions are good. Didn't seem to help at all...

* fastfit: with large data, fit on subset of data first, then go to all data

* Add documentation for kernels, esp. factor ones

* Better handling of ignore input when using factors?

* Transform inputs

* Speed up Cubic, triangle, ratquad, periodic, powerexp k/grad

* no param est gave error, dC_dparam can't be calculated

* estimate_run_time(). Run it when creating object. Give estimated time if long.

* factor trend, LM ignore inds

* transform X, Z

* penalties in optim, on trend/kernels

* Clean up trend. b_est/m_est, jitter, aren't used

* All restarts had error: doesn't say what error is, so sometimes it's
an error the user should know. Like when I test useCM.

* Knowledge gradient: multiple starts for optim

* When giving in formula/data
  * fix plot for 1D, 2D, marginal, marginal random, kernel
  * add to doc
  * message what chosen kernel is

* qEI with mixopt:
    * picks same point multiple times b/c of mean uncertainty
    * error converting mopar back and forth. mopar_converted?
    * Slow to convert to df each time. Convert mopar to mopar_converted.
    * doesn't spread out much. Better with AugEI or CorEI?
    * test, doc.
    * CL or pred?

* Corrected EI: eps, tdf

* Augmented EI: eps, tdf

* Optimize any function. Avoid reimplementing maxEI, maxAugEI, max___.
  * Can use gr or fngr
  * Can do matrix eval
  * mopar or mopar converted

* Large Z variance is bad. Extend range of s2. Warning to rescale. Tell kernel
that initial optim values should be large.

* Add deprecated for old model: lifecycle::

* Check kernel start par outside of bounds set by user. Test if it gives all errors.

* Noisy EI?

* Make gpkm doc look good. Look into R7.

* mixopt needs to be able to eval all at once in multistart

* sparse pseudo-input GP. See Snelson 2006. Inherit. Change pred_one_matrix,
deviance, etc.

* Student-t process: Inherit, change pred and deviance.
