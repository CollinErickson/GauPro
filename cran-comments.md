I received an email from Brian Ripley on 2025-11-13 to fix the issues on
https://cran.r-project.org/web/checks/check_results_GauPro.html before
2025-11-27.
It was one test that failed on multiple systems.
I made changes to tests to avoid the error.

I submitted a version on 2025-11-18 that was rejected due to a documentation
error. I have fixed it.


## Test environments
* local Windows 11 install, R 4.5.2
* R-hub builder (multiple)
* Ubuntu via GitHub Actions
* Win-builder (devel and release)
* Mac-builder

## R CMD check results

(Note to self: check Rhub with rhub::rhub_check(), then 1,3,5,20,22,30)

* local Windows 11 (11/19/25): 0 errors/warnings/notes

* local Windows 11, _R_CHECK_DEPENDS_ONLY_=TRUE (11/19/25): 0 errors/warnings/notes

* GitHub Actions, Ubuntu (11/19/25): OK

* R-Hub
  intel (11/19/25): OK
  mkl (11/19/25): OK
  linux (R-devel) (11/19/25): OK
  macos (R-devel) (11/19/25): OK
  ubuntu-release (11/19/25): 1 NOTE for slow test
  windows (R-devel) (11/19/25): OK

* Win-Builder, devel (11/19/25): OK

* Win-Builder, release (11/19/25): OK

* macOS builder (11/19/25): OK

## Downstream dependencies

* comparer: This is another one of my packages. I checked it on my
laptop and it was OK. The only code change to this package was a test, so it
shouldn't affect anything else.
