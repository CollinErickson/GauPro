Minor updates to the package.
Updating this so I can get my other package, IGP,
back onto CRAN.

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

## Downstream dependencies

None.
