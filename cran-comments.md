My update was accepted to cran on 3/28/21, but then shortly after that I
received an email saying that I need to fix the errors on 
https://cran.r-project.org/web/checks/check_results_GauPro.html.

The error on there is on r-patched-solaris-x86:
  error: call of overloaded ‘log(int)’ is ambiguous

I had Rcpp code with log(10). It isn't sure if it should return an int
or float/double. I changed the code to be log(10.0), which some answers
on Stack Overflow said would fix the error.

There was also a warning on r-devel-windows-x86_64-gcc10-UCRT:
    checking if this is a source package ... WARNING
    Subdirectory 'GauPro/src' contains apparent object files/libraries
      GauPro_init.o
    Object files/libraries should not be included in a source package.

I don't fully understand this. There is no file GauPro_init.o in /src that I can see.
If there is

## Test environments
* local Windows install, R 4.0.3
* R-hub builder
* win-builder
* local Ubuntu 20.04.2 LTS, R 4.0.3
* Ubuntu 16.04.6 LTS, R 4.0.2 (Travis)

## R CMD check results

There were no ERRORs on R 4.0.3 on Windows.
There is 1 warning for qpdf not being installed, but I wasn't able to get
qpdf to work in that laptop. Since I didn't have this warning on the others,
I think it's fine.
There is 1 note.

"
File 'GauPro/libs/x64/GauPro.dll':
  Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'

It is good practice to register native routines and to disable symbol
search.
"

This message always seems to be there,
so I don't think it's a problem.

On Travis there is one note

"
installed size is 11.0Mb
  sub-directories of 1Mb or more:
    libs   9.5Mb
"

I think this is a Unix problem because I don't get it on Windows,
and I've had this problem for previous packages.

The only note on win-builder is 
"checking CRAN incoming feasibility- Maintainer."

On my personal Ubuntu 20.04.2, there are 0 errors and 0 warnings, but 1 note

"
checking installed package size ... NOTE
    installed size is 14.1Mb
    sub-directories of 1Mb or more:
      libs  12.6Mb
"

On R-hub builder Windows Server 2008 R2 SP1, R-devel, 32/64 bit, there is 1 NOTE.
"
  * checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Collin Erickson <collinberickson@gmail.com>'
  Version contains large components (0.2.3)
"

On R-hub builder	Ubuntu Linux 20.04.1 LTS, R-release, GCC, there are three notes,
but none of them are problematic.

"
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Collin Erickson <collinberickson@gmail.com>’

Version contains large components (0.2.3)
* checking installed package size ... NOTE
  installed size is 14.1Mb
  sub-directories of 1Mb or more:
    libs  12.6Mb
* checking examples ... NOTE
Examples with CPU (user + system) or elapsed time > 5s
                         user system elapsed
GauPro_kernel_model_LOO 2.896  0.021   8.736
"

On win-builder, it says:
  Status: OK

## Downstream dependencies

None.
