## Comments

* New Note: `installed size is 5.5Mb`. 
  * This is due to datasets in this package. I believe they are essential to learning all the features
  of this package.
* Added import to package `progress`. It allows better information to the user.

## Test environments
* Rhub
  * Windows Server 2008 R2 SP1, R-devel
  * Fedora Linux, R-devel, clang, gfortran
  * Ubuntu Linux 16.04 LTS, R-release, GCC
* Travis-CI
  * Ubuntu Linux 14.04.5 LTS, R-oldrel, R-release, R-devel, GCC
  * Mac OS X 10.13.3, R-oldrel, R-release, xcode 9.4.1
* Win-builder
  * R-release, R-devel

## R CMD check results

`0 errors | 0 warnings | 1 note`

## Downstream dependencies

* No reverse dependencies yet

## Known Issues (a.k.a NOTES)

* Uses the superseded package: `doSNOW`
  * `doSNOW` has a property that allows to use progress bar that `parallel` does not.
  * Working in finding a better solution to drop this dependency. Not found yet.
  
* (possibly) invalid URLs: https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html
  * Debian: libcurl throws an error on certificate check. Nothing to do about this.

* Authors@R field gives persons with non-standard roles
  * These non-standard roles where appropriately chosen using [MARC Code List for Relators](https://www.loc.gov/marc/relators/relaterm.html)
