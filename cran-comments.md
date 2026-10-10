# CRAN Comments

## Comments for this version (0.4.16)

C++11 flags are deprecated and defaults to C++17, so this requirement was removed from the package.
Also, some old links were fixed.

The new NOTE "‘-Wdate-time’ ‘-Werror=format-security’ ‘-Wformat’" appears only in my devcontainer,
not on tested environments below, nevertheless seems ok

## Test environments

- rhub:
  -  3 [VM] macos          R-* (any version)                     macos-13 on GitHub
  -  9 [CT] clang-ubsan    R-devel (2025-11-30 r89082)           Ubuntu 22.04.5 LTS
  - 10 [CT] clang16        R-devel (2025-11-29 r89077)           Ubuntu 22.04.5 LTS
  - 16 [CT] gcc-asan       R-devel (2025-11-30 r89082)           Fedora Linux 40 (Container Image)
  - 27 [CT] ubuntu-clang   R-devel (2025-11-30 r89082)           Ubuntu 22.04.5 LTS
  - 30 [CT] ubuntu-release R-4.5.2 (2025-10-31)                  Ubuntu 24.04.3 LTS
  - 31 [CT] valgrind       R-devel (2025-11-30 r89082)           Fedora Linux 38 (Container Image)

- win-builder: devel, release, oldrel

## R CMD check results

── R CMD check results ──────────────── tsmp 0.4.16 ────
Duration: 2m 3.1s

❯ checking installed package size ... NOTE
    installed size is  9.7Mb
    sub-directories of 1Mb or more:
      data   4.6Mb
      libs   3.7Mb

❯ checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

❯ checking compilation flags used ... NOTE
  Compilation used the following non-portable flag(s):
    ‘-Wdate-time’ ‘-Werror=format-security’ ‘-Wformat’

0 errors ✔ | 0 warnings ✔ | 3 notes ✖

## Downstream dependencies

* No reverse dependencies yet

## Known Issues (a.k.a NOTES)

* Found the following (possibly) invalid file URI:
  URI: .github/CODE_OF_CONDUCT.md
  From: README.md
  * This is ok.

* GNU make is a SystemRequirements.
  * Requirement of package RcppParallel.  I haven't find a workaround to solve this NOTE.

* Installed size is X Mb.
  * This is due to datasets in this package.  I believe they are essential to learning all the features of this package.

* Uses the superseded package: `doSNOW`
  * `doSNOW` has a property that allows to use progress bar that `parallel` does not.
  * Working in finding a better solution to drop this dependency.  Not found yet.

* (possibly) invalid URLs: https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html
  * Debian: libcurl throws an error on certificate check.  Nothing to do about this.

* Authors@R field gives persons with non-standard roles
  * These non-standard roles where appropriately chosen using
    [MARC Code List for Relators](https://www.loc.gov/marc/relators/relaterm.html)

## Old comments

## Comments for this version (0.4.15)

- [x] Fixed CRAN issues; CRAN packages with clang-UBSAN errors (mail from Prof.  Brian Ripley on 1st August):
