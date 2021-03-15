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
There were no ERRORs R 4.0.3.
There is 1 warning for qpdf not being installed, but I wasn't able to get
qpdf to work in that laptop.
There is 1 note.

"
File 'GauPro/libs/x64/GauPro.dll':
  Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'

It is good practice to register native routines and to disable symbol
search.
"

I have tried everything to get rid of this note,
but can't get it to go away.
I have done what I think is the correct thing
based on all the comments I've read online.
Also I don't get this message on Travis or Win-builder,
so I don't think it is a real problem.

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

On my Ubuntu 20.04.2, there are 0 errors and 0 warnings, but 1 note

"
checking installed package size ... NOTE
    installed size is 14.1Mb
    sub-directories of 1Mb or more:
      libs  12.6Mb
"

## Downstream dependencies

None.
