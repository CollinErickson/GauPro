This package was reaccepted to CRAN on 2/12/23. I received an email from Brian
Ripley on 2/15/23 telling me to correct the error on 
https://cran.r-project.org/web/checks/check_results_GauPro.html.

The error on that page was from a unreliable test. I have improved the test to
make it much more reliable.

## Test environments
* local Windows 11 install, R 4.2.2
* R-hub builder
* Ubuntu 20.04.5 via GitHub Actions
* Win-builder
* Mac-builder

## R CMD check results

* local Windows 11 (2/20/23): no errors/warnings/notes

* Mac-builder (2/11/23): 1 NOTE for "sub-directories of 1Mb or more", but it is expected.

* GitHub  (2/11/23): 1 NOTE for large sub-directories, but no
warnings/errors.

* R-Hub Windows Server 2022 (2/20/23): 2 NOTEs, none are real issues

* R-hub Ubuntu Linux 20.04.1 (2/9/23): NOTE, none are real issues

* R-Hub Fedora Linux (2/9/23): NOTE for sub-directory size. NOTE for a slow example.
NOTE for the HTML version of the manual. NOTE for spelling, but it's correct.

* R-Hub Debian Linux (2/9/23): OK

* Win-Builder, devel (2/20/23): OK

* Win-Builder, release (2/11/23): NOTE, not a real problem

## Downstream dependencies

None.
