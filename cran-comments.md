I got an email from Brian Ripley on Jan 11, 2023 saying that I needed to fix
the issues online by Jan 25, 2023. I think I have fixed the issues because
I don't see any problems on the following checks.


## Test environments
* local Windows 11 install, R 4.2.2
* R-hub builder
* Ubuntu 20.04.5 via GitHub Actions
* Win-builder
* Mac-builder

## R CMD check results

* local Windows 11 (2/11/23): no errors/warnings/notes

* Mac-builder (2/11/23): 1 NOTE for "sub-directories of 1Mb or more", but it is expected.

* GitHub  (2/11/23): 1 NOTE for large sub-directories, but no
warnings/errors.

* R-Hub Windows Server 2022 (2/9/23): NOTE, none are real issues

* R-hub Ubuntu Linux 20.04.1 (2/9/23): NOTE, none are real issues

* R-Hub Fedora Linux (2/9/23): NOTE for sub-directory size. NOTE for a slow example.
NOTE for the HTML version of the manual. NOTE for spelling, but it's correct.

* R-Hub Debian Linux (2/9/23): OK

* Win-Builder, devel (2/11/23): NOTE, not a real problem.

* Win-Builder, release (2/11/23): NOTE, not a real problem

## Downstream dependencies

None.
