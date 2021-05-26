# CRAN Comments

## Comments

Fixed CRAN issues:

-   Please replace `\dontrun{}` by `\donttest{}` in your Rd-files.
-   Using `tempdir()` in test functions
-   Added `on.exit(par(old_par))` and similar

## Test environments

-   GitHub Actions (ubuntu-16.04): devel, release, oldrel
-   GitHub Actions (windows): release, oldrel
-   GitHub Actions (macos): release
-   rhub: all platforms (*nix, solaris, macos, windows, with gcc and clang), all ok,
  except:
    -   `solaris-x86-patched-ods` in which a dependency fails, not this package.
    -   `linux-x86_64-centos-epel`: not available: 'spelling'; found 'abort'; found 'printf'. I could not find why the code behaves like this on this distribution.
    -   `ubuntu-rchk`: Bioconductor does not yet build and check packages for R version 4.2; see <https://bioconductor.org/install>. There is nothing in the package code that asks for R 4.2 version.
-   Raspbian (armv7l-linux-gnueabihf, 32-bits): release
-   win-builder: devel, release, oldrel

## R CMD check results

`0 errors | 0 warnings | 4 notes`

## Downstream dependencies

-   No reverse dependencies yet

## Current comments

-   None yet.

## Known Issues (a.k.a NOTES)

-   [x] Found the following (possibly) invalid file URI:
  URI: .github/CODE_OF_CONDUCT.md From: README.md

-   [x] GNU make is a SystemRequirements.
    -   Requirement of package RcppParallel. I haven't find a workaround to solve this NOTE.

-   [x] Installed size is 7.5Mb.
    -   This is due to datasets in this package. I believe they are essential to learning all the features of this package.

-   [x] Uses the superseded package: `doSNOW`
    -   `doSNOW` has a property that allows to use progress bar that `parallel` does not.
    -   Working in finding a better solution to drop this dependency. Not found yet.

-   [x] (possibly) invalid URLs: <https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html>
    -   Debian: libcurl throws an error on certificate check. Nothing to do about this.

-   [x] Authors@R field gives persons with non-standard roles
    -   These non-standard roles where appropriately chosen using [MARC Code List for Relators](https://www.loc.gov/marc/relators/relaterm.html)
