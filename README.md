README
================
Francisco Bischoff
\- 23 Oct 2018

<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/logo.png" align="right" style="float:right;" />

# Time Series with Matrix Profile

[![Packagist](https://img.shields.io/badge/license-GPL--3-brightgreen.svg)](https://choosealicense.com/licenses/gpl-3.0/)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![CRAN
version](http://www.r-pkg.org/badges/version/tsmp)](https://cran.r-project.org/package=tsmp)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/tsmp)](https://cran.r-project.org/package=tsmp)

|               | Build                                                                                                                                                                        | Dev                                                                                                                                                                           |
| ------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Linux x86\_64 | [![Build Status](https://travis-ci.com/franzbischoff/tsmp.svg?branch=master)](https://travis-ci.com/franzbischoff/tsmp)                                                      | [![Build Status](https://travis-ci.com/franzbischoff/tsmp.svg?branch=develop)](https://travis-ci.com/franzbischoff/tsmp)                                                      |
| OSX           | [![Build Status](https://travis-ci.com/franzbischoff/tsmp.svg?branch=master)](https://travis-ci.com/franzbischoff/tsmp)                                                      | [![Build Status](https://travis-ci.com/franzbischoff/tsmp.svg?branch=develop)](https://travis-ci.com/franzbischoff/tsmp)                                                      |
| Windows       | [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/franzbischoff/tsmp?branch=master&svg=true)](https://ci.appveyor.com/project/franzbischoff/tsmp) | [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/franzbischoff/tsmp?branch=develop&svg=true)](https://ci.appveyor.com/project/franzbischoff/tsmp) |
| Coverage      | [![codecov](https://codecov.io/gh/franzbischoff/tsmp/branch/master/graph/badge.svg)](https://codecov.io/gh/franzbischoff/tsmp)                                               | [![codecov](https://codecov.io/gh/franzbischoff/tsmp/branch/develop/graph/badge.svg)](https://codecov.io/gh/franzbischoff/tsmp)                                               |

## Overview

R Functions implementing UCR Matrix Profile Algorithm
(<http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>).

This package allows you to use the Matrix Profile concept as a toolkit.

This package provides:

  - Algorithms to build a Matrix Profile: STAMP, STOMP, SCRIMP++, SIMPLE
    and MSTOMP.
  - Algorithms for MOTIF search for Unidimensional and Multidimensional
    Matrix Profiles.
  - Algorithm for Chains search for Unidimensional Matrix Profile.
  - Algorithms for Semantic Segmentation (FLUSS) and Weakly Labeled data
    (SDTS).
  - Algorithm for Salient Subsections detection allowing MDS plotting.
  - Basic plotting for all outputs generated here.
  - Sequencial workflow, see below.

<!-- end list -->

``` r
# Basic workflow:
matrix <- tsmp(data, window_size = 30) %>% find_motif(n_motifs = 3) %>% plot()

# SDTS still have a unique way to work:
model <- sdts_train(data, labels, windows)
result <- sdts_predict(model, data, round(mean(windows)))
```

Please refer to the [User
Manual](https://franzbischoff.github.io/tsmp/reference/) for more
details.

Please be welcome to suggest
improvements.

### Performance on an Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz using a random walk dataset

``` r
set.seed(2018)
data <- cumsum(sample(c(-1, 1), 40000, TRUE))
```

|               | Elapsed Time | Data size | Window size | Threads |
| ------------- | :----------: | :-------: | :---------: | :-----: |
| `stomp_par()` |    52.72s    |   40000   |    1000     |    8    |
| `scrimp()`    |    92.44s    |   40000   |    1000     |    1    |
| `stomp()`     |   133.16s    |   40000   |    1000     |    1    |
| `stamp_par()` |   140.25s    |   40000   |    1000     |    8    |
| `stamp()`     |   262.03s    |   40000   |    1000     |    1    |

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
  - SCRIMP (single-thread, not for AB-joins yet)
  - Time Series Chains
  - Multivariate STOMP (mSTOMP)
  - Multivariate MOTIF Search (from mSTOMP)
  - Salient Subsequences search for Multidimensional Space
  - Scalable Dictionary learning for Time Series (SDTS) prediction
  - FLUSS (Fast Low-cost Unipotent Semantic Segmentation)
  - SiMPle-Fast (Fast Similarity Matrix Profile for Music Analysis and
    Exploration)
  - Annotation vectors (e.g., Stop-word MOTIF bias, Actionability bias)
  - FLUSS Arc Plot and SiMPle Arc Plot
  - Misc:
      - MASS v2.0
      - Fast moving average
      - Fast moving SD

## Roadmap

  - Exact Detection of Variable Length Motifs
  - Profile-Based Shapelet Discovery
  - GPU-STOMP
  - Real-time version of previous algorithms (STAMPI, FLOSS, etc.)
  - MASS Extensions (ADP, WQ, QwG)

## Other projects with Matrix Profile

  - Python: <https://github.com/ZiyaoWei/pyMatrixProfile>
  - Python: <https://github.com/jbeleno/owlpy>
  - Python: <https://github.com/javidlakha/matrix-profile>
  - CUDA: <https://github.com/zpzim/STOMPSelfJoin>
  - CUDA: <https://github.com/zpzim/SCAMP>

## Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](CODE_OF_CONDUCT.md). By participating in this project, you
agree to abide by its terms.
