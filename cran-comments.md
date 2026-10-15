# CRAN Comments

## Comments for this version (0.5.0)

- Architectural modernization: `tsmp` is now a pure R package (`NeedsCompilation: no`).
  All core computational algorithms (MPX, STOMP, STAMP, SCRIMP, MASS, mSTOMP) are delegated to
  the dedicated high-performance engine package `matrixprofiler`.
- Eradicated superseded dependencies: `doSNOW` and `foreach` have been completely removed.
- Removed all compiled C++ code and GNU Make requirements (`RcppParallel`/TBB eliminated).
- Standardized invalid index conventions to `NA_integer_`.
- All 363 unit tests passed cleanly with 0 failures.
- `R CMD check --as-cran` passed with 0 errors, 0 warnings, and 0 notes.
- Fixed moved/broken URLs in README.md.

## Test environments

- rhub: linux-devel, windows-devel, macos-devel, clang-ubsan (failed not because of tsmp, it is a r-hub problem), valgrind

- win-builder: devel, release

## R CMD check results

─  checking installed package size ... INFO
     installed size is  5.5Mb
     sub-directories of 1Mb or more:
       data   4.6Mb

── R CMD check results ────────────tsmp 0.5.0 ────
Duration: 38.9s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Downstream dependencies

* No reverse dependencies yet

## Known Issues (a.k.a NOTES)

* Found the following (possibly) invalid file URI:
  URI: .github/CODE_OF_CONDUCT.md
  From: README.md
  * This is ok.

* GNU make is a SystemRequirements:
  * Solved! GNU make and RcppParallel removed; package is now pure R.

* Installed size is ~5.3 Mb.
  * This is due to datasets in this package. They are essential to learning all features of this package.

* Uses the superseded package: `doSNOW`:
  * Solved! `doSNOW` and `foreach` have been completely removed. Multithreading is now handled by `matrixprofiler`.

* (possibly) invalid URLs: https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html
  * Debian: libcurl throws an error on certificate check.  Nothing to do about this.

* Authors@R field gives persons with non-standard roles
  * These non-standard roles where appropriately chosen using
    [MARC Code List for Relators](https://www.loc.gov/marc/relators/relaterm.html)

## Old comments

## Comments for this version (0.4.15)

- [x] Fixed CRAN issues; CRAN packages with clang-UBSAN errors (mail from Prof.  Brian Ripley on 1st August):
