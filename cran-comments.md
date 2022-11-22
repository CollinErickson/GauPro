This package was accepted to CRAN on 11/14/22.
The next day I received an email saying I had until 11/29/22 to fix
some build errors or the package will be removed from CRAN.

The errors were from flaky tests that use randomness and fail a small portion
of the time. I have rewritten these tests so that they will pass consistently.

I also added a reference to a related paper of mine in the DESCRIPTION
and made some other improvements.


## Test environments
* local Windows 11 install, R 4.2.2
* R-hub builder
* Ubuntu 20.04.5 via GitHub Actions
* Win-builder
* Mac-builder

## R CMD check results

* local Windows 11: no errors/warnings/notes

* Mac-builder: 1 NOTE for "sub-directories of 1Mb or more", but it is expected.

* GitHub Actions: 1 NOTE for large sub-directories, but no
warnings/errors.

* R-Hub Windows Server 2022: OK

* R-hub Ubuntu Linux 20.04.1: NOTE for sub-directory size and another
NOTE for a slow example.

* R-Hub Fedora Linux: NOTE for sub-directory size. NOTE for a slow example.
NOTE for the HTML version of the manual. There was an issue building a vignette,
but I haven't seen this issue anywhere else.

* R-Hub Debian Linux: OK

* Win-Builder: NOTE for incoming feasibility.

## Downstream dependencies

None.
