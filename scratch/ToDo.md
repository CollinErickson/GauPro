# GauPro to do

* Take in data/formula input in addition to matrix

* Calculate EI: nopt, test, doc

* max_qEI: test, doc, exact solution. Change CL to CLpred

* With many points, just estimate param from subset of data, but use all data
at end for K/Kinv so predictions are good.

* Improve kernel model plot

* progress bar for restarts

* Fewer restarts in high dim

* Give kernel in as string. If missing, give better error or just use Mat52.

* Reduce L-BFGS-B tolerance? Need to change tests too.

* Add documentation for kernels, esp. factor ones

* Better handling of ignore input when using factors?

* Transform inputs

* plot2D option to do se instead of mean, or side by side

* optim NaN starting value (found on 1D doing EI)

* plot LOO calibration

* Speed up m52 grad, Cubic k/grad

* 3 kernel product. prod/sum need to get s2_est from children

* no param est gave error, dC_dparam can't be calculated

* Add test with repeated X. Add test with big nugget and make sure deriv still matches.

* Fix travis after success

* Prevent negative s2, e.g., borehole

* optim: start by checking n points, then optimizing from best?
