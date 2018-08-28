README
================
Francisco Bischoff
\- 28 Aug 2018

<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/logo.png" style="float: right;" />

# Time Series with Matrix Profile

[![Packagist](https://img.shields.io/packagist/l/doctrine/orm.svg)](https://choosealicense.com/licenses/mit)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![CRAN
version](http://www.r-pkg.org/badges/version/tsmp)](https://cran.r-project.org/package=tsmp)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/tsmp)](https://cran.r-project.org/package=tsmp)

|               | Build                                                                                                                          | Dev                                                                                                                             |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------- |
| Linux x86\_64 | [![Build Status](https://travis-ci.com/franzbischoff/tsmp.svg?branch=master)](https://travis-ci.com/franzbischoff/tsmp)        | [![Build Status](https://travis-ci.com/franzbischoff/tsmp.svg?branch=develop)](https://travis-ci.com/franzbischoff/tsmp)        |
| OSX           | [![Build Status](https://travis-ci.com/franzbischoff/tsmp.svg?branch=master)](https://travis-ci.com/franzbischoff/tsmp)        | [![Build Status](https://travis-ci.com/franzbischoff/tsmp.svg?branch=develop)](https://travis-ci.com/franzbischoff/tsmp)        |
| Coverage      | [![codecov](https://codecov.io/gh/franzbischoff/tsmp/branch/master/graph/badge.svg)](https://codecov.io/gh/franzbischoff/tsmp) | [![codecov](https://codecov.io/gh/franzbischoff/tsmp/branch/develop/graph/badge.svg)](https://codecov.io/gh/franzbischoff/tsmp) |

## Overview

R Functions implementing UCR Matrix Profile Algorithm
(<http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>).

This is under development and is intended to be a general purpose MP
toolkit

After basic tools are finished and API is mature, further functions for
‘Classification’, ‘MOTIF extraction’, ‘MDS visualization’ etc. will be
added.

Please be welcome to suggest improvements.

### Performance on an Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz

|                | Elapsed Time | Data size | Window size | Cores |
| -------------- | :----------: | :-------: | :---------: | :---: |
| `stomp.par()`  |    45.17s    |   55000   |     150     |   8   |
| `stomp()`      |    76.68s    |   55000   |     150     |   1   |
| `mstomp.par()` |   113.03s    |   55000   |     150     |   8   |
| `mstomp()`     |   238.81s    |   55000   |     150     |   1   |
| `stamp.par()`  |   852.11s    |   55000   |     150     |   8   |
| `stamp()`      |   2862.79s   |   55000   |     150     |   1   |

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
  - Time Series Chains
  - Multivariate STOMP (mSTOMP)
  - Multivariate MOTIF Search (from mSTOMP)
  - Salient Subsequences search for Multidimensional Space
  - Scalable Dictionary learning for Time Series (SDTS) prediction
  - FLUSS (Fast Low-cost Unipotent Semantic Segmentation)
  - SiMPle-Fast (Fast Similarity Matrix Profile for Music Analysis and
    Exploration)
  - Misc:
      - MASS v2.0
      - Fast moving average
      - Fast moving SD

## Road map

  - FLUSS Arc Plot and SiMPle Arc Plot
  - Annotation vectors (e.g.: Stop-word MOTIF bias, Actionability bias)
  - MOTIFs under Uniform Scaling
  - GPU-STOMP
  - Real-time version of previous algorithms (STAMPI, FLOSS, etc)
  - MASS Extensions (ADP, WQ, QwG)
  - SCRIMP (waiting for publication)

## Other projects with Matrix Profile

  - Python: <https://github.com/ZiyaoWei/pyMatrixProfile>
  - Python: <https://github.com/jbeleno/owlpy>
  - Python: <https://github.com/javidlakha/matrix-profile>

## Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree
to abide by its terms.
