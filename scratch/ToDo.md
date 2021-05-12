# GauPro to do

* Factor kernel: example

* Ordered factor kernel: example

* Ignore input kernel: example

* Take in data/formula input in addition to matrix

* Calculate EI: nopt, test, doc

* max_qEI: test, doc, exact solution

* With many points, just estimate param from subset of data, but use all data
at end for K/Kinv so predictions are good.

* Change print, summary

* Improve kernel model plot

* progress bar for restarts

* Fewer restarts in high dim

* Give kernel in as string. If missing, give better error or just use Mat52.

* Reduce L-BFGS-B tolerance? Need to change tests too.

* Add documentation for kernels, esp. factor ones

* Better handling of ignore input when using factors?

* Add cubic correlation

* Transform inputs

* plot2D option to do se instead of mean, or side by side

* optim NaN starting value (found on 1D doing EI)
