## Package update
This is an update of the TRexSelector R package:

  * Version is increased to 1.0.0 (stable version).
  * Supporting scientific papers were updated.
  * New methods and extensions of the T-Rex selector were added (e.g., dependency-aware T-Rex, Screen-T-Rex, etc.).
  * DESCRIPTION, NEWS, README were updated.
  * Package was polished (vignettes, etc.).
  

## Resubmission
This is a resubmission that addresses the comments on the initial submission. The following changes were made:

 * All unit tests use no more than two cores. In the initial submission, some unit tests used more than two cores. However, these unit tests were skipped on CRAN (testthat::skip_on_cran()). In the current submission, no unit tests are skipped on CRAN.
 * Functions trex() and random_experiments() were polished.

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs.

## Downstream dependencies
There are currently no downstream dependencies for this package.
