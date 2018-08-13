README
================
Francisco Bischoff
August 13, 2018

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Time Series Matrix-Profile <img src="man/figures/logo.png" align="right" />

Main:
[![Build Status main](https://travis-ci.com/franzbischoff/tsmp.svg?branch=master)](https://travis-ci.com/franzbischoff/tsmp)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/r-lib/tsmp?branch=master&svg=true)](https://ci.appveyor.com/project/r-lib/tsmp)
[![codecov](https://codecov.io/gh/franzbischoff/tsmp/branch/master/graph/badge.svg)](https://codecov.io/gh/franzbischoff/tsmp)
[![CRAN
version](http://www.r-pkg.org/badges/version/tsmp)](https://cran.r-project.org/package=tsmp)
Dev:
[![Build Status dev](https://travis-ci.com/franzbischoff/tsmp.svg?branch=feature%2Fcode_status)](https://travis-ci.com/franzbischoff/tsmp)
[![codecov](https://codecov.io/gh/franzbischoff/tsmp/branch/feature/code_status/graph/badge.svg)](https://codecov.io/gh/franzbischoff/tsmp)

R Functions implementing UCR Matrix Profile Algorithm
(<http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>).

This is under development and is intended to be a MP toolkit

Further files will be provided with actual use cases as
‘classification’, ‘motif extraction’, ‘MDS visualization’, etc.

Please be welcome to suggest improvements.

## Overview

Testing your code can be painful and tedious, but it greatly increases
the quality of your code. **tsmp** tries to make testing as fun as
possible, so that you get a visceral satisfaction from writing tests.
Testing should be addictive, so you do it all the time. To make that
happen, tsmp:

  - Provides functions that make it easy to describe what you expect a
    function to do, including catching errors, warnings, and messages.

  - Easily integrates in your existing workflow, whether it’s informal
    testing on the command line, building test suites, or using R CMD
    check.

  - Displays test progress visually, showing a pass, fail, or error for
    every expectation. If you’re using the terminal or a recent version
    of RStudio, it’ll even colour the output.

tsmp draws inspiration from the xUnit family of testing packages, as
well as from many of the innovative ruby testing libraries, like
[rspec](http://rspec.info/), [testy](https://github.com/ahoward/testy),
[bacon](https://github.com/chneukirchen/bacon) and
[cucumber](https://cucumber.io).

tsmp is the most popular unit testing package for R and is used by
thousands of CRAN packages.

If you’re not familiar with tsmp, the [testing
chapter](http://r-pkgs.had.co.nz/tests.html) in [R
packages](http://r-pkgs.had.co.nz/) gives a good overview, along with
workflow advice and concrete examples.

## Installation

``` r
# Install the released version from CRAN
install.packages("tsmp")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("franzbischoff/tsmp")
```

# Python implementation

<https://github.com/ZiyaoWei/pyMatrixProfile>
