# GauPro to do

* EI/CorEI/AugEI/qEI: nopt, test, doc
  * Use t-dist

* Add documentation for kernels, esp. factor ones

* Add readme/documentation for trends

* Speed up triangle, ratquad, periodic, powerexp k/grad

* Knowledge gradient: multiple starts for optim

* Large Z variance is bad. Extend range of s2. Warning to rescale. Tell kernel
that initial optim values should be large.

* Plot2D:
  add axis names, either X1/X2 or colnames (need to add option to CF::gcf_grid)
  fix for factors

* Check kernel start par outside of bounds set by user. Test if it gives all errors.

* Make gpkm doc look good. Look into R7.

* sparse pseudo-input GP. See Snelson 2006. Inherit. Change pred_one_matrix,
deviance, etc.

* Student-t process: Inherit, change pred and deviance.

* kernel should take in X, Z and set params based on that. E.g., s2 max is
diff(range(Z))^2. Or else give message to normalize.

* Standardize for X.

* factor kernel: clean up jitter/start
