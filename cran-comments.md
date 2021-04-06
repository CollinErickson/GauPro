(Update for 4/5/21)

My last submission still had two valgrind errors in the examples.
I have fixed these and am resubmitting. Sorry again.

(Update for 4/4/21)

Uwe emailed me with the details of the valgrind error upon my last resubmission.
I finally found that there was a single perilous example with bad inputs.
I fixed the example and added a test so it won't happen again.

I reran the check on my computer, Travis, R-hub, and win-devel
with no new issues.

I apologize for these failed resubmissions. I hate to think I wasted your time.
I have tried to fix all errors each time.

(From previous submission)

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
I added .o, .so, and .dll files to my .Rbuildignore, but I'm not that will fix it.
I couldn't recreate this error on any other system.
Since it's just a warning and only shows up on one system, I think it's okay to resubmit
since this will fix the known error.

I reran R CMD check on my computer, on R-hub, and win-devel and had no new issues come up.


I resubmitted this on 4/1/21, but it did not go through because I hadn't fixed the valgrind error.
I specifically checked the valgrind output before submitting it, but on the valgrind 00check.log 
(https://www.stats.ox.ac.uk/pub/bdr/memtests/valgrind/GauPro/00check.log), 
the status says OK with no notes, warnings, or errors, so I thought it was fine.

On checking the Ex.Rout file, I see that there are errors on 3 lines of code. 
It appears that all of them are when I use "exp" in Rcpp code. According to this Stack Overflow by Dirk 
(https://stackoverflow.com/questions/40997722/rcpp-armadillo-rstudio-says-exp-is-ambiguous),
I think the issue is that I just used "exp" instead of "std:exp" or "arma::exp".

I have fixed these errors and also ran on R-hub with valgrind, which had no errors,
so it is good to go. I also reran on R-hub and win-devel, with no change in results
on any of those.




## Test environments
* local Windows install, R 4.0.3
* R-hub builder
* win-builder
* local Ubuntu 20.04.2 LTS, R 4.0.3
* Ubuntu 16.04.6 LTS, R 4.0.2 (Travis)
* R-hub with valgrind

## R CMD check results

There were no ERRORs on R 4.0.3 on Windows.
There is 1 warning for qpdf not being installed, but I wasn't able to get
qpdf to work in that laptop. Since I didn't have this warning on the others,
I think it's fine.
There is 1 NOTE:
    File 'GauPro/libs/x64/GauPro.dll':
      Found no calls to: 'R_registerRoutines', 'R_useDynamicSymbols'
    
    It is good practice to register native routines and to disable symbol
    search.

This message always seems to be there,
so I don't think it's a problem.

On Travis there is one NOTE:
    installed size is 11.0Mb
      sub-directories of 1Mb or more:
        libs   9.5Mb

I think this is a Unix problem because I don't get it on Windows,
and I've had this problem for previous packages.

The only note on win-builder is 
"checking CRAN incoming feasibility- Maintainer."

On my personal Ubuntu 20.04.2, there are 0 errors and 0 warnings, but 1 note:
    checking installed package size ... NOTE
        installed size is 14.1Mb
        sub-directories of 1Mb or more:
          libs  12.6Mb

On R-hub builder Fedora Linux, R-devel, clang, gfortran, there are 2 NOTEs:
    * checking CRAN incoming feasibility ... NOTE
    Maintainer: ‘Collin Erickson <collinberickson@gmail.com>’
    
    Days since last update: 4
    * checking installed package size ... NOTE
      installed size is  7.4Mb
      sub-directories of 1Mb or more:
        libs   6.0Mb

On R-hub builder Windows Server 2008 R2 SP1, R-devel, 32/64 bit, there is 1 NOTE.
    * checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Collin Erickson <collinberickson@gmail.com>'
    
    Days since last update: 4

On R-hub builder Ubuntu Linux 20.04.1 LTS, R-release, GCC, there are two notes,
but neither of them are problematic.
    * checking CRAN incoming feasibility ... NOTE
    Maintainer: ‘Collin Erickson <collinberickson@gmail.com>’
    
    Days since last update: 4
    * checking installed package size ... NOTE
      installed size is 14.1Mb
      sub-directories of 1Mb or more:
        libs  12.6Mb

On win-builder, it says:
    * checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Collin Erickson <collinberickson@gmail.com>'
    
    Days since last update: 4

Because of the valgrind error, I also ran with valgrind on R-hub using rhub::check_with_valgrind(),
where there was only a standard note.

    checking installed package size ... NOTE
    
    installed size is 12.8Mb
    
    sub-directories of 1Mb or more:
    
    libs 11.4Mb

## Downstream dependencies

None.
