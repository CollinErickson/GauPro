I made some updates to the package.
I fixed a test that gave CRAN issues recently.


## Test environments
* local Windows 11 install, R 4.4.3
* R-hub builder (multiple)
* Ubuntu 24.04.2 via GitHub Actions
* Win-builder (devel and release)
* Mac-builder

## R CMD check results

* local Windows 11 (3/9/25): no errors/warnings/notes

* local Windows 11, _R_CHECK_DEPENDS_ONLY_=TRUE (3/9/25): no errors/warnings/notes

* Mac-builder (3/9/25): 1 NOTE
1 NOTE for "sub-directories of 1Mb or more", but it is expected.

* GitHub Actions, Ubuntu (3/9/25): 1 NOTE for large sub-directories, but no
warnings/errors.

* R-Hub ubuntu gcc-12 (3/9/25): 1 NOTE for slow tests

* R-Hub linux (R-devel) (3/9/25): OK

* R-Hub macos (R-devel) (3/9/25): OK

* R-Hub windows (R-devel) (3/9/25): OK

* R-Hub intel (3/9/25): OK

* Win-Builder, devel (3/9/25): OK

* Win-Builder, release (3/9/25): OK

## Downstream dependencies

* (formerly) comparer: This is another one of my packages. I checked it on my
laptop and it was OK.
