# GauPro to do

* EI: nopt, test, doc
  * Use t-dist
  * Use dEI_dx in maxEI

* Add documentation for kernels, esp. factor ones

* Better handling of ignore input when using factors?

* Speed up triangle, ratquad, periodic, powerexp k/grad

* Knowledge gradient: multiple starts for optim

* When giving in formula/data
  * add to doc
  * message what chosen kernel is

* qEI with mixopt:
    * picks same point multiple times b/c of mean uncertainty
    * error converting mopar back and forth. mopar_converted?
    * Slow to convert to df each time. Convert mopar to mopar_converted.
    * doesn't spread out much. Better with AugEI or CorEI?
    * test, doc.
    * CL or pred?

* Corrected EI: tdf

* Augmented EI: tdf

* Optimize any function. Avoid reimplementing maxEI, maxAugEI, max___.
  * Can use gr or fngr
  * Can do matrix eval
  * mopar or mopar converted

* Large Z variance is bad. Extend range of s2. Warning to rescale. Tell kernel
that initial optim values should be large.

* Plot2D:
  add axis names, either X1/X2 or colnames
  fix for factors

* Add deprecated for old model: lifecycle::

* Check kernel start par outside of bounds set by user. Test if it gives all errors.

* Make gpkm doc look good. Look into R7.

* mixopt needs to be able to eval all at once in multistart

* sparse pseudo-input GP. See Snelson 2006. Inherit. Change pred_one_matrix,
deviance, etc.

* Student-t process: Inherit, change pred and deviance.

* kernel should take in X, Z and set params based on that. E.g., s2 max is
diff(range(Z))^2. Or else give message to normalize.

* Check maxqEI for different EI type, see if different

* plotmarginal w/ factors needs level names

* GowerFactorKernel: 1 if equal, s_i if not. Param for each dim
* HammingFactorKernel: 1 if equal, s if not, single param?

* pseudo r-sq should use df

* factor kernel: clean up logp, offdiagequal, jitter/start
