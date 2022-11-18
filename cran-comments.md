This package was accepted to CRAN on 11/14/22.
The next day I received an email saying I had until 11/29/22 to fix
some build errors or the package will be removed from CRAN.

The errors are from flaky tests that use randomness and fail a small portion
of the time. I have rewritten these tests so that they will pass consistently.



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

* R-hub Ubuntu: NOTE for sub-directory size, NOTE for slow examples,
NOTE for new submission.

## Downstream dependencies

None.
