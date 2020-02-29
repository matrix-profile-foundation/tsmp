README
================
Francisco Bischoff
\- 28 Feb 2020

<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/logo.png" align="right" style="float:right;" />

# Time Series with Matrix Profile

<!-- badges: start -->
[![Packagist](https://img.shields.io/badge/license-GPL--3-brightgreen.svg)](https://choosealicense.com/licenses/gpl-3.0/)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![CRAN
version](http://www.r-pkg.org/badges/version/tsmp)](https://cran.r-project.org/package=tsmp)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/tsmp)](https://cran.r-project.org/package=tsmp)
[![CircleCI build
status](https://circleci.com/gh/matrix-profile-foundation/tsmp.svg?style=svg)](https://circleci.com/gh/matrix-profile-foundation/tsmp)
<!-- badges: end -->

|               | Build                                                                                                                                                                             | Dev                                                                                                                                                                                 |
| ------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Linux x86\_64 | [![Build Status](https://travis-ci.com/matrix-profile-foundation/tsmp.svg?branch=master)](https://travis-ci.com/matrix-profile-foundation/tsmp)                                   | [![Build Status](https://travis-ci.com/matrix-profile-foundation/tsmp.svg?branch=develop)](https://travis-ci.com/matrix-profile-foundation/tsmp)                                    |
| OSX           | [![Build Status](https://travis-ci.com/matrix-profile-foundation/tsmp.svg?branch=master)](https://travis-ci.com/matrix-profile-foundation/tsmp)                                   | [![Build Status](https://travis-ci.com/matrix-profile-foundation/tsmp.svg?branch=develop)](https://travis-ci.com/matrix-profile-foundation/tsmp)                                    |
| Windows       | [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/byfyqncr60ten98g/branch/master?svg=true)](https://ci.appveyor.com/project/franzbischoff/tsmp/branch/master) | [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/byfyqncr60ten98g/branch/develop?svg=true)](https://ci.appveyor.com/project/franzbischoff/tsmp/branch/develop) |
| Coverage      | [![codecov](https://codecov.io/gh/matrix-profile-foundation/tsmp/branch/master/graph/badge.svg)](https://codecov.io/gh/matrix-profile-foundation/tsmp)                            | [![codecov](https://codecov.io/gh/matrix-profile-foundation/tsmp/branch/develop/graph/badge.svg)](https://codecov.io/gh/matrix-profile-foundation/tsmp)                             |

## Overview

R Functions implementing UCR Matrix Profile Algorithm
(<http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>).

This package allows you to use the Matrix Profile concept as a toolkit.

This package provides:

  - Algorithms to build a Matrix Profile: STAMP, STOMP, SCRIMP++,
    SIMPLE, MSTOMP and VALMOD.
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
matrix <- tsmp(data, window_size = 30) %>%
  find_motif(n_motifs = 3) %T>%
  plot()

# SDTS still have a unique way to work:
model <- sdts_train(data, labels, windows)
result <- sdts_predict(model, data, round(mean(windows)))
```

Please refer to the [User
Manual](https://matrix-profile-foundation.github.io/tsmp/reference/) for
more details.

Please be welcome to suggest improvements.

### Performance on an Intel(R) Core(TM) i7-7700 CPU @ 3.60GHz using a random walk dataset

``` r
set.seed(2018)
data <- cumsum(sample(c(-1, 1), 40000, TRUE))
```

#### Current version benchmark

|             | Elapsed Time(s) | Data Size | Window Size | Threads | Lang |
| ----------- | --------------: | --------: | ----------: | ------: | :--- |
| `mpx_par`   |            0.59 |     40000 |        1000 |       8 | Rcpp |
| `mpx`       |            1.94 |     40000 |        1000 |       1 | Rcpp |
| `stomp_par` |           38.90 |     40000 |        1000 |       8 | R    |
| `stomp`     |           85.13 |     40000 |        1000 |       1 | R    |
| `scrimp`    |          123.07 |     40000 |        1000 |       1 | R    |
| `stamp_par` |          925.45 |     40000 |        1000 |       8 | R    |
| `stamp`     |         3776.86 |     40000 |        1000 |       1 | R    |

## Installation

``` r
# Install the released version from CRAN
install.packages("tsmp")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("matrix-profile-foundation/tsmp")
```

## Currently available Features

  - STAMP (single and multi-thread versions)
  - STOMP (single and multi-thread versions)
  - STOMPi (On-line version)
  - SCRIMP (single-thread, not for AB-joins yet)
  - Time Series Chains
  - Multivariate STOMP (mSTOMP)
  - Multivariate MOTIF Search (from mSTOMP)
  - Salient Subsequences search for Multidimensional Space
  - Scalable Dictionary learning for Time Series (SDTS) prediction
  - FLUSS (Fast Low-cost Unipotent Semantic Segmentation)
  - FLOSS (Fast Low-cost On-line Unipotent Semantic Segmentation)
  - SiMPle-Fast (Fast Similarity Matrix Profile for Music Analysis and
    Exploration)
  - Annotation vectors (e.g., Stop-word MOTIF bias, Actionability bias)
  - FLUSS Arc Plot and SiMPle Arc Plot
  - Exact Detection of Variable Length Motifs (VALMOD)
  - MPdist: Matrix Profile Distance
  - Time Series Snippets
  - Subsetting Matrix Profiles (`head()`, `tail()`, `[`, etc.)
  - Misc:
      - MASS v2.0
      - MASS v3.0
      - MASS extensions: ADP (Approximate Distance Profile, with PAA)
      - MASS extensions: WQ (Weighted Query)
      - MASS extensions: QwG (Query with Gap)
      - Fast moving average
      - Fast moving SD

## Roadmap

  - Profile-Based Shapelet Discovery
  - GPU-STOMP

## Other projects with Matrix Profile

  - Python: <https://github.com/target/matrixprofile-ts>
  - Python: <https://github.com/ZiyaoWei/pyMatrixProfile>
  - Python: <https://github.com/jbeleno/owlpy>
  - Python: <https://github.com/javidlakha/matrix-profile>
  - Python: <https://github.com/shapelets/khiva-python>
  - R: <https://github.com/shapelets/khiva-r>
  - Matlab: <https://github.com/shapelets/khiva-matlab>
  - Java: <https://github.com/shapelets/khiva-java>
  - Java: <https://github.com/ensozos/Matrix-Profile>
  - Kotlin: <https://github.com/shapelets/khiva-kotlin>
  - C++ (CUDA and OPENCL): <https://github.com/shapelets/khiva>
  - CUDA: <https://github.com/zpzim/STOMPSelfJoin>
  - CUDA: <https://github.com/zpzim/SCAMP>

## Matrix Profile Foundation

Our next step unifying the Matrix Profile implementation in several
programming languages.

Visit: [Matrix Profile Foundation](https://matrixprofile.org)

## Code of Conduct

Please note that the ‘tsmp’ project is released with a [Contributor Code
of
Conduct](https://github.com/matrix-profile-foundation/tsmp/blob/master/.github/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.
