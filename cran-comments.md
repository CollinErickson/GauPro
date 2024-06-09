I fixed the warning on CRAN.
I received an email that I had to fix these issues by 5/28 or it would be
removed from CRAN, and then another email last saying that I had until 6/20.

I'm not sure I'm using the R-Hub checks the right way now that there are 24
choices. I picked a handful on a variety of platforms, I hope that is good enough.

Mac builder (`devtools::check_mac_release()`) didn't work for me, but I ran both
of the MacOS options in R-Hub with no issues. The Mac build link just gives
a 404 error.

## Test environments
* local Windows 11 install, R 4.4.0
* R-hub builder
* Ubuntu 22.04.4 via GitHub Actions
* Win-builder
* (didn't work this time: Mac-builder)

## R CMD check results

* local Windows 11 (6/2/24): no errors/warnings/notes

* Mac-builder (4/10/23): 1 NOTE for "sub-directories of 1Mb or more", but it is expected.

* GitHub (6/2/24): 1 NOTE for large sub-directories, but no
warnings/errors.

* R-Hub linux (ubuntu-latest on GitHub) (6/7/24): 1 NOTE for sub-directory size.
with my package.

* R-Hub atlas (Fedora Linux 38) (6/7/24): OK

* R-Hub clang19 (Ubuntu 22.04.4 LTS) (6/9/24): 1 NOTE for
  Compilation used the following non-portable flag(s):
    ‘-Wp,-D_FORTIFY_SOURCE=3’

* R-Hub gcc14 (Fedora Linux 40) (6/9/24): OK

* R-Hub ubuntu release (Ubuntu 22.04.4 LTS) (6/9/24): 1 NOTE for
"Running R code in ‘testthat.R’ had CPU time 3 times elapsed time"

* R-Hub macOS (macos-13 on GitHub) (6/9/24): OK

* R-Hub macOS-arm64 (macos-latest on GitHub) (6/9/24): OK

* Win-Builder, devel (6/3/24): 1 NOTE for a possibly invalid URL, but the URL
is fine.

* Win-Builder, release (6/4/24): 1 NOTE for a possibly invalid URL, but the URL
is fine.

## Downstream dependencies

* comparer: This is another one of my packages. The changes in this package
won't affect it since I only changed the tests. I also ran R CMD CHECK on
on comparer and it passed with no issues.
