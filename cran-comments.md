I received an email from Brian
Ripley on 3/29/23 telling me to correct the error on 
https://cran.r-project.org/web/checks/check_results_GauPro.html.

The error on that page was from a unreliable test. I have improved the test to
make it reliable.

## Test environments
* local Windows 11 install, R 4.2.3
* R-hub builder
* Ubuntu 22.04.2 via GitHub Actions
* Win-builder
* Mac-builder

## R CMD check results

* local Windows 11 (4/7/23): no errors/warnings/notes

* Mac-builder (4/7/23): 1 NOTE for "sub-directories of 1Mb or more", but it is expected.

* GitHub (4/8/23): 1 NOTE for large sub-directories, but no
warnings/errors.

* R-Hub Windows Server 2022 (failed 4/4/23 and 4/8): It had a problem loading
standard packages, it's not the fault of my package.

* R-Hub Fedora Linux (4/8/23): NOTE for sub-directory size. NOTE for a slow example.

* R-Hub Debian Linux (4/8/23): PREPERROR

* R-Hub Ubuntu Linux (4/8/23): NOTEs for sub-directory size and slow example.

* Win-Builder, devel (4/8/23): OK

* Win-Builder, release (4/8/23): NOTE for a slow example

## Downstream dependencies

* comparer: This is another one of my packages. The changes in this package
won't affect it since I only changed the tests. I also ran R CMD CHECK on
on comparer and it passed with no issues.
