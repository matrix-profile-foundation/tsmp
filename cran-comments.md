## Comments
Fixed CRAN issues:
- CRAN packages with clang-UBSAN errors (mail from Prof.  Brian Ripley on 1st August)

## Test environments
* Rhub
  * Windows Server 2022, R-devel, 64 bit - all ok
  * Ubuntu Linux 20.04.1 LTS, R-release, GCC - all ok
  * Ubuntu Linux 20.04.1 LTS, R-devel, GCC - all ok
  * Fedora Linux, R-devel, GCC - all ok
  * Fedora Linux, R-devel, clang, gfortran - all ok
  * Debian Linux, R-release, GCC ASAN/UBSAN - all ok
  * Debian Linux, R-patched, GCC - all ok
  * Debian Linux, R-devel, GCC, no long double - all ok
  * Debian Linux, R-devel, GCC ASAN/UBSAN - all ok
  * Debian Linux, R-devel, clang, ISO-8859-15 locale - all ok
  * Apple Silicon (M1), macOS 11.6 Big Sur, R-release - all ok
  * macOS 10.13.6 High Sierra, R-release, brew - all ok
  * Oracle Solaris 10, x86, 32 bit, R release, Oracle Developer Studio 12.6 - all ok (except for `testthat` not
    available)
* Win-builder
  * R-release, R-oldrelease, R-devel - all ok

## R CMD check results

`0 errors | 0 warnings | 2 notes`

* GNU make is a SystemRequirements.
  * Requirement of package RcppParallel.  I haven't find a workaround to solve this NOTE.

* Installed size is X Mb.
  * This is due to datasets in this package.  I believe they are essential to learning all the features of this package.

## Downstream dependencies

* No reverse dependencies yet

## Known Issues (a.k.a NOTES)

* Found the following (possibly) invalid file URI:
  URI: .github/CODE_OF_CONDUCT.md
  From: README.md
  * This is ok.

* GNU make is a SystemRequirements.
  * Requirement of package RcppParallel.  I haven't find a workaround to solve this NOTE.

* Installed size is 7.5Mb.
  * This is due to datasets in this package.  I believe they are essential to learning all the features of this package.

* Uses the superseded package: `doSNOW`
  * `doSNOW` has a property that allows to use progress bar that `parallel` does not.
  * Working in finding a better solution to drop this dependency.  Not found yet.

* (possibly) invalid URLs: https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html
  * Debian: libcurl throws an error on certificate check.  Nothing to do about this.

* Authors@R field gives persons with non-standard roles
  * These non-standard roles where appropriately chosen using
    [MARC Code List for Relators](https://www.loc.gov/marc/relators/relaterm.html)

