This package was removed from CRAN on 10/3/22 because I didn't update it.
Sorry for taking so long to get back to this.



## Test environments
* local Windows 11 install, R 4.2.2
* R-hub builder
* Ubuntu 20.04.5 via GitHub Actions
* Win-builder
* Mac-builder

## R CMD check results

* local Windows 11: no errors/warnings/notes

* Mac-builder: 1 NOTE for "sub-directories of 1Mb or more", but it is expected.
Also an error for a failed test, but it's a test with randomness that passes
on other attempts.

* GitHub Actions: Also the 1 NOTE for large sub-directories, but no
warnings/errors.

## Downstream dependencies

None.
