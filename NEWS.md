# tsmp (development version)

NEWS
================
Francisco Bischoff
\- 28 Feb 2020

<!-- NEWS.md is generated from NEWS.Rmd. Please edit that file -->

# tsmp 0.4.8

  - Added MPdist algorithm to compare two time series.
  - Added `find_snippet()` that uses MPdist to show representative data.
  - Added `mpx()` algorithm that doesn’t depends on FFT.
  - Added `pmp()` pan-matrix profile.
  - Improvement with several implementations in Rcpp. Computation speed
    is much faster.
  - Added `compute()`, `analyze()` and `visualize()`. Starting point in
    the unified API from MPF.

# tsmp 0.3.5

  - Changed `mass()` to `dist_profile()`, including options to different
    algorithms.
  - Added MASS\_V3 and MASS\_Weighted to `dist_profile()` algorithms.
  - Function `dist_profile()` allows Query with Gap (QwG) and
    Approximate Distance Profile (ADP, with PAA)
  - Fixed long runtime of FFT for some data sizes, using MASS\_V3. Issue
    \#36.
  - Added `stompi_update()` that updates the matrix profile allowing
    real-time computation.
  - Added `floss()` which can do real-time FLUSS computation.
  - Added subset operator `[` for `tsmp` objects.
  - Added `tail()` and `head()` for `tsmp` objects.
  - Changed `find_motif()` for `MultiMatrixProfile` to report Motifs
    correctly.

# tsmp 0.3.4

  - Added `find_discord()` and its `print()` and `plot()` functions.
  - Changed `plot()` for motifs to show where are the neighbors. Same
    for discord.
  - Added `valmod()` for Variable Length Motif Discovery.
  - Changed `find_motif()` for compatibility with `valmod()`.

# tsmp 0.3.3

  - Fixed `find_chains()` not returning the longest chain. Issue \#33

# tsmp 0.3.2

  - Fixed Matrix Profile print, dimensions are now reported correctly.
  - Fixed pipe imports. Issue \#22
  - Fixed bug with `vars`. Issue \#23
  - Changed package license to GPL-3.
  - Changed verbose mode, added one more step to separate messages from
    progression bar.
  - Fixed SCRIMP and added PRE-SCRIMP, so this is the SCRIMP++. AB-join
    not yet implemented.
  - Changed progress bar for a better one from `progress` package.
  - Added Print and Plot to SiMPle. Issue \#24
  - Added Print and Plot to Salient.

# tsmp 0.3.1

## IMPORTANT

  - This version is a complete restructuration. The API has changed, and
    the workflow is more friendly. This API is intended to be stable,
    and from now on any change will pass through the “Deprecated” stage.

## Added Features

  - Outputs have a prettier print format.
  - Outputs have a plot function. Try to plot a `tsmp()` output for
    example.
  - Now functions can work in `%>%` (pipe), e.g. `tsmp() %>%
    find_motif()`. Except for SDTS that has a proper way to work.
  - Added a wrapper function called `tsmp()` that handles the several
    algorithms available.
  - Added `as.*` functions to allow you to switch classes if you want,
    e.g.: `as.matrixprofile()`.
  - Changed all functions from dotted.case to snake\_case (except `as.*`
    functions).
  - Added Annotation Vectors.
  - Fixed STOMP crash with Joins.
  - Added support to query \< data in Joins.
  - SCRIMP (experimental).

# tsmp 0.2.15

  - Code linting.
  - Added Salient Subsequences search.

# tsmp 0.2.14

  - Added SiMPle (Fast Similarity Matrix Profile for Music Analysis and
    Exploration).
  - Added FLUSS (Fast Low-cost Unipotent Semantic Segmentation).
  - Added \[find\_chains()\] to look for chains primitives.
  - Added Multivariate MOTIF Search (from mSTOMP)
  - Changed dependency from beepr to audio (actually beepr depends on
    audio, so fewer dependencies).
  - Added a `NEWS.md` file to track changes to the package.

# tsmp 0.2.12

  - Added Multivariate STOMP parallel version.
  - Added SDTS algorithm (Scalable Dictionary learning for Time Series).

# tsmp 0.1.0

  - STAMP and STAMP parallel Algorithm.
  - Multivariate STOMP algorithm.
  - MASS algorithm.
