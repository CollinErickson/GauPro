# GauPro 0.2.17.9000

# GauPro 0.2.17

Fixed a test that caused error on CRAN.

Changed the class name of GauPro_base() from "GauPro"" to "GauPro_base" to avoid
documentation error.

Accepted to CRAN on 11/21/25.

# GauPro 0.2.16

Fixed a test that caused error on CRAN.

Accepted to CRAN on 8/26/25.

# GauPro 0.2.15

Fixed a test that caused error on CRAN tests.

Accepted to CRAN on 4/8/25.

# GauPro 0.2.14

Bug fix from predictions when there are no categorical predictors and formula
input is used. 

Added isotropic option for Gaussian, Exponential, Matern 3/2, Matern 5/2,
and Triangle kernels.

Combined README and main vignette for better documentation.

Accepted to CRAN on 3/9/25.

# GauPro 0.2.13

Fixed Suggests issue to get back on CRAN. Packages in Suggests were moved
to Depends, removed, or protected by requireNamespace.

Accepted to CRAN on 9/26/24.

# GauPro 0.2.12

Fixed CRAN warning.

Added `k_xyz(...)` alias for kernels (replaces `xyz$new(...)`)

Accepted to CRAN on 6/10/24.

# GauPro 0.2.11

Fixed unreliable test to keep it on CRAN (again).

# GauPro 0.2.8

Fixed unreliable test to keep it on CRAN.

Accepted to CRAN on 2/27/23.

# GauPro 0.2.7

Improved summary, importance, plots.

Added gradpredvar, AugmentedEI, CorrectedEI, optimize_fn.

GauPro was removed from CRAN on 1/25/23, this puts it back on CRAN.
Accepted to CRAN on 2/12/23.

# GauPro 0.2.6

GP kernel model maxEI can now be run using mixopt to account for discrete
inputs.

Improved GP kernel model workability when input is formula and data frame.

Added importance for kernel model, greatly improved summary.

Fixed error in tests to prevent being removed from CRAN.

Accepted to CRAN on 11/24/22.

# GauPro 0.2.5

Added kernels for factors.

Changed default number of restarts to zero, and added checking more starting
points. Should make it faster.

Can give in data as data and formula instead of matrix and vector.

Package was removed from CRAN on 10/3/22. This fixes the issue.

Accepted to CRAN on 11/14/22.

# GauPro 0.2.4

Very minor change in Rcpp code to remove CRAN error.

Accepted to CRAN on 4/11/21.

# GauPro 0.2.3

Minor changes.

Accepted to CRAN on 3/28/21, but was notified of error on same day.

# GauPro 0.2.2

Fixing Valgrind error from 0.2.1.

Accepted to CRAN on 9/11/2017.

# GauPro 0.2.1

Fixing minor errors from the 0.2.0 version.


# GauPro 0.2.0

Added kernel models that use kernels and trends.

# GauPro 0.1.0

Releasing for the first time.

Accepted by CRAN on 10/11/16
