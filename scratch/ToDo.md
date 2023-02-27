# GauPro to do

* EI/CorEI/AugEI/qEI: nopt, test, doc
  * Use t-dist

* Add documentation for kernels, esp. factor ones

* Speed up triangle, ratquad, periodic, powerexp k/grad

* Knowledge gradient: multiple starts for optim

* When giving in formula/data
  * add to doc
  * message what chosen kernel is

* Large Z variance is bad. Extend range of s2. Warning to rescale. Tell kernel
that initial optim values should be large.

* Plot2D:
  add axis names, either X1/X2 or colnames
  fix for factors

* Add deprecated for old model: lifecycle::

* Check kernel start par outside of bounds set by user. Test if it gives all errors.

* Make gpkm doc look good. Look into R7.

* sparse pseudo-input GP. See Snelson 2006. Inherit. Change pred_one_matrix,
deviance, etc.

* Student-t process: Inherit, change pred and deviance.

* kernel should take in X, Z and set params based on that. E.g., s2 max is
diff(range(Z))^2. Or else give message to normalize.

* Standardize for X.

* plotmarginal w/ factors needs level names

* factor kernel: clean up jitter/start

* Fix GHA/Linux bug with cubic

* increasing nugget to get invertibility, only print last

* badger badges
