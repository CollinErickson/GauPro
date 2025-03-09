I made some updates to the package and I fixed a test that gave CRAN issues
recently.


## Test environments
* local Windows 11 install, R 4.4.3
* R-hub builder (multiple)
* Ubuntu 22.04.5 via GitHub Actions
* Win-builder (devel and release)
* Mac-builder

## R CMD check results

* local Windows 11 (3/9/25): no errors/warnings/notes

* local Windows 11, _R_CHECK_DEPENDS_ONLY_=TRUE (3/9/25): no errors/warnings/notes

* Mac-builder (9/22/24): 1 NOTE, 1 WARNING
1 NOTE for "sub-directories of 1Mb or more", but it is expected.
1 WARNING for package ‘mixopt’ was built under R version 4.4.1. Not a real issue.

* GitHub Actions, Ubuntu (9/22/24): 1 NOTE for large sub-directories, but no
warnings/errors.

* R-Hub ubuntu gcc-12 (9/22/24): 1 NOTE for a slow test

* R-Hub linux (R-devel) (9/22/24): 1 NOTE for a large sub-directory

* R-Hub macos (R-devel) (9/22/24): OK

* R-Hub windows (R-devel) (9/22/24): OK

* Win-Builder, devel (9/22/24): 1 NOTE for "checking CRAN incoming feasibility".
New submission, possibly misspelled words that are fine, possibly broken link
that is fine.

* Win-Builder, release (9/22/24):  1 NOTE for "checking CRAN incoming feasibility".
New submission, possibly misspelled words that are fine, possibly broken link
that is fine.

## Downstream dependencies

* (formerly) comparer: This is another one of my packages, but it was also removed from 
CRAN for the Suggests issue. I'll fix this next.
