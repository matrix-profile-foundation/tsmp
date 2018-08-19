README
================
Francisco Bischoff
\- 19 Aug 2018

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Time Series with Matrix Profile <img src="man/figures/logo.png" align="right" />

[![Packagist](https://img.shields.io/packagist/l/doctrine/orm.svg)](https://choosealicense.com/licenses/mit)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Build
Status](https://travis-ci.com/franzbischoff/tsmp.svg?branch=develop)](https://travis-ci.com/franzbischoff/tsmp)
[![codecov](https://codecov.io/gh/franzbischoff/tsmp/branch/develop/graph/badge.svg)](https://codecov.io/gh/franzbischoff/tsmp)
[![CRAN
version](http://www.r-pkg.org/badges/version/tsmp)](https://cran.r-project.org/package=tsmp)

## Overview

R Functions implementing UCR Matrix Profile Algorithm
(<http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>).

This is under development and is intended to be a general purpose MP
toolkit

After basic tools are finished and API is mature, further functions for
‘Classification’, ‘MOTIF extraction’, ‘MDS visualization’ etc. will be
added.

Please be welcome to suggest improvements.

## Installation

``` r
# Install the released version from CRAN
install.packages("tsmp")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("franzbischoff/tsmp")
```

## Currently available Features

  - STAMP (single and multi-thread versions)
  - STOMP (single and multi-thread versions)
  - Multivariate STOMP (mSTOMP)
  - Multivariate MOTIF Search (from mSTOMP)
  - Scalable Dictionary learning for Time Series (SDTS) prediction
  - FLUSS (Fast Low-cost Unipotent Semantic Segmentation)
  - Misc:
      - MASS v2.0
      - Fast moving average
      - Fast moving SD

## Road map

  - Improve Joins outputs
  - Multidimensional Space
  - Time Series Chains
  - FLUSS Arc Plot and SiMPle Arc Plot
  - Annotation vectors (e.g.: Stop-word MOTIF bias, Actionability bias)
  - SiMPle-Fast (Fast Similarity Matrix Profile for Music Analysis and
    Exploration)
  - MOTIFs under Uniform Scaling
  - GPU-STOMP
  - Real-time version of previous algorithms (STAMPI, FLOSS, etc)
  - MASS Extensions (ADP, WQ, QwG)
  - SCRIMP (waiting for publication)

## Other projects with Matrix Profile

  - Python: <https://github.com/ZiyaoWei/pyMatrixProfile>

## Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree
to abide by its terms.
