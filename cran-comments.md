I received an email from Brian
Ripley on 3/29/23 telling me to correct the error on 
https://cran.r-project.org/web/checks/check_results_GauPro.html.

The error on that page was from a unreliable test. I have improved the test to
make it reliable.

## Test environments
* local Windows 11 install, R 4.2.3
* R-hub builder
* Ubuntu 20.04.5 via GitHub Actions
* Win-builder
* Mac-builder

## R CMD check results

* local Windows 11 (4/7/23): no errors/warnings/notes

* Mac-builder (4/7/23): 1 NOTE for "sub-directories of 1Mb or more", but it is expected.

* GitHub  (2/25/23): 1 NOTE for large sub-directories, but no
warnings/errors.

* R-Hub Windows Server 2022 (2/26/23, failed 4/4): NOTE, none are real issues

* R-Hub Fedora Linux (4/4/23): NOTE for sub-directory size. NOTE for a slow example.

* R-Hub Debian Linux (2/21/23, prep error 4/4): OK

* Win-Builder, devel (2/20/23): OK

* Win-Builder, release (2/25/23): NOTE, not a real problem

## Downstream dependencies

None.
