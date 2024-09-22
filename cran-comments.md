I received an email from Brian Ripley on 8/27 that I needed to fix this package.
One of my other packages, ContourFunctions, had been removed from CRAN,
and this package was failing noSuggests. I have resubmitted ContourFunctions.
I fixed this package so that everything in Suggests was moved to Depends
or is only used within requireNamespace.


## Test environments
* local Windows 11 install, R 4.4.0
* R-hub builder
* Ubuntu 22.04.5 via GitHub Actions
* Win-builder
* (didn't work this time: Mac-builder)

## R CMD check results

* local Windows 11 (9/21/24): no errors/warnings/notes

* local Windows 11, _R_CHECK_DEPENDS_ONLY_=TRUE (9/22/24): no errors/warnings/notes

* Mac-builder (9/22/24): 1 NOTE, 1 WARNING
1 NOTE for "sub-directories of 1Mb or more", but it is expected.
1 WARNING for package ‘mixopt’ was built under R version 4.4.1. Not a real issue.

* GitHub Actions, Ubuntu (9/21/24): 1 NOTE for large sub-directories, but no
warnings/errors.

* R-Hub linux (ubuntu-latest on GitHub) (6/7/24): 1 WARNING, 1 NOTE
1 NOTE for for sub-directory size.
Warning: package ‘mixopt’ was built under R version 4.4.1

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
