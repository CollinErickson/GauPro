# GauPro to do

* Replace main vignette GauPro() to gpkm(), follow GH doc

* Fix preds negative: mean up to 0, preds up to nug*s2hat

* Fix error from CRAN

* EI/CorEI/AugEI/qEI: nopt, test, doc
  * Use t-dist

* Add documentation for kernels, esp. factor ones

* Add readme/documentation for trends

* Speed up triangle, ratquad, periodic, powerexp k/grad

* Knowledge gradient: multiple starts for optim

* Plot2D:
  add axis names, either X1/X2 or colnames (need to add option to CF::gcf_grid)
  fix for factors

* Make gpkm doc look good. Look into R7.

* sparse pseudo-input GP. See Snelson 2006. Inherit. Change pred_one_matrix,
deviance, etc.

* Student-t process: Inherit, change pred and deviance.

* kernel should take in X, Z and set params based on that. E.g., s2 max is
diff(range(Z))^2. Or else give message to normalize.

* Standardize for X.

* change kernels to functions k_: documentation, tests, examples, readme

* change trends from R6 to function t_

* gpk$plot1D(ymax=1000) doesn't use ymax?
