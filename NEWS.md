NEWS
================
Francisco Bischoff
\- 01 Oct 2018

<!-- NEWS.md is generated from NEWS.Rmd. Please edit that file -->

# tsmp 0.3.2

  - Fixed Matrix Profile print, dimensions are now reported correctly.
  - Fixed pipe imports. Issue \#22
  - Fixed bug with `vars`. Issue \#23
  - Changed package license to GPL-3.
  - Changed verbose mode, added one more step to separate messages from
    progression bar.

# tsmp 0.3.1

## IMPORTANT

  - This version is a complete reestructuration. API has changed and
    workflow is more friendly. This API is intended to be stable and
    from now one any change will pass throught the “Deprecated” stage.

## Added Features

  - Outputs have a prettier print format.
  - Outputs have a plot function. Try to plot a `tsmp()` output for
    example.
  - Now functions can work in `%>%` (pipe), e.g. `tsmp() %>%
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
    audio, so less dependencies).
  - Added a `NEWS.md` file to track changes to the package.

# tsmp 0.2.12

  - Added Multivariate STOMP parallel version.
  - Added SDTS algorithm (Scalable Dictionary learning for Time Series).

# tsmp 0.1.0

  - STAMP and STAMP parallel Algorithm.
  - Multivariate STOMP algorithm.
  - MASS algorithm.
