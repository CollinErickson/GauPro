Version 0.2.1 was accepted 9/8/17 but there was a
valgrind error. I believe I have found and fixed
the issue.

## Test environments
* local Windows install, R 3.4.1
* win-builder
* Ubuntu 12.04.5 LTS, R 3.4.1 (Travis)

## R CMD check results
There were no ERRORs or WARNINGs on Windows using R 3.4.1.
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
installed size is  6.6Mb
  sub-directories of 1Mb or more:
    libs   6.1Mb
"

I think this is a Unix problem because I don't get it on Windows,
and I've had this problem for previous packages.

The only note on win-builder is 
"checking CRAN incoming feasibility- Maintainer."

## Downstream dependencies

None.
