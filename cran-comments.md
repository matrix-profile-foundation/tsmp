## Comments

* Changed dependency from beepr to audio (actually beepr depends on audio, so less dependencies).
* Added persons to DESCRIPTION as their code/research was implemented in this package.
* Changed dotted.case to snake_case as usual in tidyverse. I've started with Google R style, but I
concluded this is better. Since the package is starting, I don't think will have any impact on users.

## Test environments
* Rhub
  * Windows Server 2008 R2 SP1, R-devel
  * Fedora Linux, R-devel, clang, gfortran
  * Ubuntu Linux 16.04 LTS, R-release, GCC
* Travis-CI
  * Ubuntu Linux 14.04.5 LTS, R-release, GCC
* Win-builder
  * R-devel
  * R-release

## R CMD check results

`0 errors | 0 warnings | 0 notes`

## Downstream dependencies

* No reverse dependencies yet

## Known Issues (a.k.a NOTES)

* Uses the superseded package: `doSNOW`
  * `doSNOW` has a property that allow to use progress bar that `parallel` does not.
  * Working in finding a good solution to drop this dependency.

* Authors@R field gives persons with non-standard roles
  * These non-standard roles where chosen properly using [MARC Code List for Relators](https://www.loc.gov/marc/relators/relaterm.html)
