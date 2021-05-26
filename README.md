README
================
Francisco Bischoff

26 mai 2021

<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/logo.png" align="right" style="float:right;"/>

## Time Series with Matrix Profile

<!-- badges: start -->

[![License](https://img.shields.io/badge/License-GPL--3.0-green.svg)](https://choosealicense.com/licenses/gpl-3.0/)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html)

[![Lint](https://github.com/matrix-profile-foundation/tsmp/workflows/Lint/badge.svg?branch=main)](https://github.com/jimhester/lintr)
[![R-CMD-check](https://github.com/matrix-profile-foundation/tsmp/workflows/R-CMD-check/badge.svg?branch=main)](https://r-pkgs.org/r-cmd-check.html)

[![codecov](https://codecov.io/gh/matrix-profile-foundation/tsmp/branch/main/graph/badge.svg?token=w7AmbwhNvn)](https://codecov.io/gh/matrix-profile-foundation/tsmp)

[![CRAN
version](https://www.r-pkg.org/badges/version/tsmp)](https://cran.r-project.org/package=tsmp)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/tsmp)](https://cran.r-project.org/package=tsmp)

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/matrix-profile-foundation/tsmp/master)

<!-- badges: end -->

## Important News!!!

The `tsmp` package is being modified to allow a great change in speed
using the `matrixprofiler` package. We will make all efforts to keep it
back compatible.

A slightly more explained text is available
[here](https://franzbischoff.rbind.io/posts/presenting-matrixprofiler-a-fast-matrix-profile-implementation-in-r/).

### Overview

R Functions implementing UCR Matrix Profile Algorithm
(<http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>).

This package allows you to use the Matrix Profile concept as a toolkit.

This package provides:

-   Algorithms to build a Matrix Profile: STAMP, STOMP, SCRIMP++,
    SIMPLE, MSTOMP and VALMOD.
-   Algorithms for MOTIF search for Unidimensional and Multidimensional
    Matrix Profiles.
-   Algorithm for Chains search for Unidimensional Matrix Profile.
-   Algorithms for Semantic Segmentation (FLUSS) and Weakly Labeled data
    (SDTS).
-   Algorithm for Salient Subsections detection allowing MDS plotting.
-   Basic plotting for all outputs generated here.
-   Sequencial workflow, see below.

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

### Installation

``` r
# Install the released version from CRAN
install.packages("tsmp")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("matrix-profile-foundation/tsmp")
```

### Currently available Features

-   STAMP (single and multi-thread versions)

-   STOMP (single and multi-thread versions)

-   STOMPi (On-line version)

-   SCRIMP (single-thread, not for AB-joins yet)

-   Time Series Chains

-   Multivariate STOMP (mSTOMP)

-   Multivariate MOTIF Search (from mSTOMP)

-   Salient Subsequences search for Multidimensional Space

-   Scalable Dictionary learning for Time Series (SDTS) prediction

-   FLUSS (Fast Low-cost Unipotent Semantic Segmentation)

-   FLOSS (Fast Low-cost On-line Unipotent Semantic Segmentation)

-   SiMPle-Fast (Fast Similarity Matrix Profile for Music Analysis and
    Exploration)

-   Annotation vectors (e.g., Stop-word MOTIF bias, Actionability bias)

-   FLUSS Arc Plot and SiMPle Arc Plot

-   Exact Detection of Variable Length Motifs (VALMOD)

-   MPdist: Matrix Profile Distance

-   Time Series Snippets

-   Subsetting Matrix Profiles (`head()`, `tail()`, `[`, etc.)

-   Misc:

    -   MASS v2.0
    -   MASS v3.0
    -   MASS extensions: ADP (Approximate Distance Profile, with PAA)
    -   MASS extensions: WQ (Weighted Query)
    -   MASS extensions: QwG (Query with Gap)
    -   Fast moving average
    -   Fast moving SD

### Roadmap

-   Profile-Based Shapelet Discovery
-   GPU-STOMP

### Other projects with Matrix Profile

-   Python: <https://github.com/target/matrixprofile-ts>
-   Python: <https://github.com/ZiyaoWei/pyMatrixProfile>
-   Python: <https://github.com/jbeleno/owlpy>
-   Python: <https://github.com/javidlakha/matrix-profile>
-   Python: <https://github.com/shapelets/khiva-python>
-   R: <https://github.com/shapelets/khiva-r>
-   Matlab: <https://github.com/shapelets/khiva-matlab>
-   Java: <https://github.com/shapelets/khiva-java>
-   Java: <https://github.com/ensozos/Matrix-Profile>
-   Kotlin: <https://github.com/shapelets/khiva-kotlin>
-   C++ (CUDA and OPENCL): <https://github.com/shapelets/khiva>
-   CUDA: <https://github.com/zpzim/STOMPSelfJoin>
-   CUDA: <https://github.com/zpzim/SCAMP>

### Matrix Profile Foundation

Our next step unifying the Matrix Profile implementation in several
programming languages.

Visit: [Matrix Profile Foundation](https://matrixprofile.org)

### Benchmarks

New benchmarks for the new package
[matrixprofiler](https://CRAN.R-project.org/package=matrixprofiler) is
available at [RPubs](https://rpubs.com/franzbischoff/matrixprofiler).

### Package dependencies

<center>

![](man/figures/dependency_plot-1.png)<!-- -->

### Code of Conduct

Please note that the tsmp project is released with a [Contributor Code
of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
