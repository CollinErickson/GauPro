I submitted this on 10/11 and received an email from Kurt Hornik
saying that it wouldn't install. I don't get this error on my 
computers, but I think I fixed it by adding a makevars file.
I was previously failing on Travis (couldn't figure out why),
but after this change it is now working.

I submitted again on 10/11 and received an email from Uwe Ligges
requesting that I single quote software names in the description. 
I have done this and am resubmitting.

## Test environments
* local Windows install, R 3.3.1
* UNIX on a cluster, R 3.1.2
* win-builder
* Ubuntu 12.04.5 LTS, R 3.3.1 (Travis)

## R CMD check results
There were no ERRORs or WARNINGs or NOTEs on Windows using R 3.3.1.

I get two notes on UNIX that don't seem like real problems (NEWS.md and installed package size).

The only note on win-builder is Maintainer- New submission.

## Downstream dependencies

None, this is the first release.
