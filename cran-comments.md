I received an email from Brian Ripley to fix the issues on
https://cran.r-project.org/web/checks/check_results_GauPro.html before 8/25/25.
I made changes to tests to avoid the error and other future potential errors
from flakey tests.


## Test environments
* local Windows 11 install, R 4.5.1
* R-hub builder (multiple)
* Ubuntu 24.04.2 via GitHub Actions
* Win-builder (devel and release)
* Mac-builder

## R CMD check results

(Note to self: check Rhub with rhub::rhub_check(), then 1,3,5,20,21)

* local Windows 11 (8/24/25): 0 errors/warnings/notes

* local Windows 11, _R_CHECK_DEPENDS_ONLY_=TRUE (8/24/25): 0 errors/warnings/notes

* Mac-builder (4/6/25): 1 NOTE
1 NOTE for "sub-directories of 1Mb or more", but it is expected.

* GitHub Actions, Ubuntu (4/6/25): 1 NOTE for large sub-directories, but no
warnings/errors.

* R-Hub
  mkl (4/6/25): OK
  linux (R-devel) (4/6/25): OK
  macos (R-devel) (4/6/25): OK
  windows (R-devel) (4/6/25): OK
  intel (4/6/25): OK

* Win-Builder, devel (8/24/25): OK

* Win-Builder, release (4/6/25): OK

## Downstream dependencies

* comparer: This is another one of my packages. I checked it on my
laptop and it was OK. The only code change to this package was a test, so it
shouldn't affect anything else.
