README
================
Francisco Bischoff
– 15 Aug 2018

# Time Series Matrix-Profile <img src="man/figures/logo.png" align="right" />

[![Build
Status](https://travis-ci.com/franzbischoff/tsmp.svg?branch=develop)](https://travis-ci.com/franzbischoff/tsmp)
[![codecov](https://codecov.io/gh/franzbischoff/tsmp/branch/develop/graph/badge.svg)](https://codecov.io/gh/franzbischoff/tsmp)

## Overview

R Functions implementing UCR Matrix Profile Algorithm
(<http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>).

This is under development and is intended to be a MP toolkit

Further files will be provided with actual use cases as
‘Classification’, ‘MOTIF extraction’, ‘MDS visualization’, etc.

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
  - Annotation vectors (e.g.: Stop-word MOTIF bias, Actionability bias)
  - SiMPle-Fast (Fast Similarity Matrix-Profile for Music Analysis and
    Exploration)
  - MOTIFs under Uniform Scaling
  - GPU-STOMP
  - Real-time version of previous algorithms (STAMPI, FLOSS, etc)
  - MASS Extensions (ADP, WQ, QwG)
  - SCRIMP (waiting for publication)

## Other projects with Matrix-Profile

<https://github.com/ZiyaoWei/pyMatrixProfile>
