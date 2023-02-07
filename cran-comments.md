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
