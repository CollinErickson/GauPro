I fixed the warning on CRAN.

## Test environments
* local Windows 11 install, R 4.2.3
* R-hub builder
* Ubuntu 22.04.2 via GitHub Actions
* Win-builder
* Mac-builder

## R CMD check results

* local Windows 11 (4/9/23): no errors/warnings/notes

* Mac-builder (4/10/23): 1 NOTE for "sub-directories of 1Mb or more", but it is expected.

* GitHub (4/9/23): 1 NOTE for large sub-directories, but no
warnings/errors.

* R-Hub Windows Server 2022 (4/9/23): NOTE, it doesn't seem to be a problem
with my package.

* R-Hub Fedora Linux (4/10/23): NOTEs for sub-directory size and slow example.

* R-Hub Debian Linux (4/10/23): PREPERROR. It says:
"ERROR: dependency ‘RcppArmadillo’ is not available for package ‘GauPro’".
This isn't the fault of my package.

* R-Hub Ubuntu Linux (4/10/23): NOTEs for sub-directory size and slow example.

* Win-Builder, devel (4/10/23): OK

* Win-Builder, release (4/8/23): NOTE for a slow example

## Downstream dependencies

* comparer: This is another one of my packages. The changes in this package
won't affect it since I only changed the tests. I also ran R CMD CHECK on
on comparer and it passed with no issues.
