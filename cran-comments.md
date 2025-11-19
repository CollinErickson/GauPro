I received an email from Brian Ripley on 2025-11-13 to fix the issues on
https://cran.r-project.org/web/checks/check_results_GauPro.html before
2025-11-17.
It was one test that failed on multiple systems.
I made changes to tests to avoid the error.


## Test environments
* local Windows 11 install, R 4.5.1
* R-hub builder (multiple)
* Ubuntu via GitHub Actions
* Win-builder (devel and release)

## R CMD check results

(Note to self: check Rhub with rhub::rhub_check(), then 1,3,5,20,21,29)

* local Windows 11 (8/24/25): 0 errors/warnings/notes

* local Windows 11, _R_CHECK_DEPENDS_ONLY_=TRUE (8/24/25): 0 errors/warnings/notes

* GitHub Actions, Ubuntu (8/25/25): OK

* R-Hub
  intel (8/25/25): OK
  mkl (8/25/25): OK
  linux (R-devel) (8/25/25): OK
  macos (R-devel) (8/25/25): OK
  ubuntu-release (8/25/25): OK

* Win-Builder, devel (8/24/25): OK

* Win-Builder, release (8/24/25): OK

## Downstream dependencies

* comparer: This is another one of my packages. I checked it on my
laptop and it was OK. The only code change to this package was a test, so it
shouldn't affect anything else.
