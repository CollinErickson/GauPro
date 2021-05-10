# GauPro to do

* Factor kernel: example

* Ordered factor kernel: example

* Ignore input kernel: example

* Take in data/formula input in addition to matrix

* Calculate EI: nopt, test, doc

* EI should use mu instead of Z for noise. Should it be sd of mean, not pred?
See http://krasserm.github.io/2018/03/21/bayesian-optimization/.

* max_qEI: test, doc, exact solution

* With many points, just estimate param from subset of data, but use all data
at end for K/Kinv so predictions are good.

* Change print, summary

* Improve kernel model plot

* Predict mean vs value?

* progress bar for restarts

* Fewer restarts in high dim

* Give kernel in as string. If missing, give better error or just use Mat52.

* Reduce L-BFGS-B tolerance? Need to change tests too.

* Add documentation for kernels, esp. factor ones

* Better handling of ignore input when using factors?

* Add cubic correlation

* Transform inputs

* Change S3 plot to R6

* plot1D show interval for mean as well as predictive

* plot2D option to do se instead of mean, or side by side
