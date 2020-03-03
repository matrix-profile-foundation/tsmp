## Comments
quick-fix: checkmate package was in suggests but is imports.
This version added Rcpp implementations and RcppParallel to allow multi-threading.
Progressivelly I'll convert the bottlenecks to Rcpp and hope to get rid of `doSNOW` for parallelization.

## Test environments
* Rhub
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  * Ubuntu Linux 16.04 LTS, R-release, GCC
  * Fedora Linux, R-devel, clang, gfortran
  * Debian Linux, R-devel, GCC ASAN/UBSAN
* Travis-CI
  * Ubuntu Linux 14.04.5 LTS, R-oldrel, R-release, R-devel, GCC
  * Mac OS X 10.13.3, R-oldrel, R-release, xcode 9.4.1
* Win-builder
  * R-release, R-devel

## R CMD check results

`0 errors | 0 warnings | 4 notes`

## Downstream dependencies

* No reverse dependencies yet

## Known Issues (a.k.a NOTES)

* Found the following (possibly) invalid file URI:
  URI: .github/CODE_OF_CONDUCT.md
  From: README.md
  * This is ok.

* GNU make is a SystemRequirements.
  * Requirement of package RcppParallel. I haven't find a workaround to solve this NOTE.

* Installed size is 7.5Mb. 
  * This is due to datasets in this package. I believe they are essential to learning all the features
    of this package.

* Uses the superseded package: `doSNOW`
  * `doSNOW` has a property that allows to use progress bar that `parallel` does not.
  * Working in finding a better solution to drop this dependency. Not found yet.
  
* (possibly) invalid URLs: https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html
  * Debian: libcurl throws an error on certificate check. Nothing to do about this.

* Authors@R field gives persons with non-standard roles
  * These non-standard roles where appropriately chosen using [MARC Code List for Relators](https://www.loc.gov/marc/relators/relaterm.html)

