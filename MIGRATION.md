# API, Compatibility, and Migration Guide (Breaking Changes)

**Packages:** `matrixprofiler` (0.1.10 $\to$ 0.1.12) & `tsmp` (0.4.16 $\to$ 0.5.0)
**Ecosystem:** Matrix Profile Foundation (MPF)
**Date:** October 2026

---

## 1. Overview of the New Architecture

In this release cycle, the Matrix Profile R ecosystem underwent a major architectural modernization:

1. **`matrixprofiler` (0.1.12)**: Serves as the **low-level computational engine** in C++17. It uses `RcppThread` for safe native parallelism (free from fork conflicts on POSIX and with no dependency on Intel TBB) and zero-allocation data structures, such as linear $O(N)$ monotonic deques and $O(K)$ algorithms.
2. **`tsmp` (0.5.0)**: Becomes a **pure R package (`NeedsCompilation: no`)**. It removes C++ compilers, Makefiles, and legacy cluster libraries (`doSNOW`, `foreach`), delegating all computationally intensive operations to `matrixprofiler`.

This clear separation delivers speedups of up to hundreds of times, straightforward installation on any operating system, and eliminates memory leaks or hangs in multicore environments.

---

## 2. `matrixprofiler`: What Changed (0.1.10 $\to$ 0.1.12)

### 2.1 New Functions Exported in the API
* **`mstomp()`**: Native C++ implementation of the Multivariate STOMP algorithm, accepting `data`, `window_size`, `must_dim`, `exc_dim`, `n_workers`, and `progress`.
* **`mpdist()`**: Native C++ implementation of MPdist (scalar `"simple"` and vector `"vector"` versions), with native multithreading support through `n_workers`.

### 2.2 Signature Changes and New Parameters
* **Sliding window functions (`windowfunc`)**:
  `mov_mean()`, `mov_sum()`, `mov_var()`, `mov_std()`, `mov_min()`, `mov_max()`, `movmean_std()`, and `muinvn()`.
  * **New parameters**: `na_strategy = c("window_na", "remove", "error")` and `na.rm = FALSE`.
* **`mass()`**:  
  * The `version` argument has been expanded to formally support: `c("v3", "v2", "weighted", "non_normalized", "absolute")`.
  * Native thread support through `n_workers` is available in all MASS variants.

---

### 2.3 Potential Breaking Changes in `matrixprofiler`

| Item | Previous Behavior (0.1.10) | New Behavior (0.1.12) | Risk & Recommended Action |
| :--- | :--- | :--- | :--- |
| **Invalid Index Convention** | In some C++ routines, windows with no match or trivial matches returned `-1L` or `0`. | **`NA_integer_`**: All index vectors and matrices now use R's canonical convention (`NA_integer_`). | ⚠️ **Medium**: If your code checked `idx == -1` or `idx <= 0`, update it to use `is.na(idx)`. |
| **`NA` Values in Moving Windows** | `mov_*` functions propagated undefined values (`NaN`), produced erratic asymptotic behavior, or overflowed. | **`na_strategy = "window_na"` (default)**: Any window that overlaps a non-finite value strictly returns `NA_real_`, without corrupting adjacent valid windows. | ⚠️ **Low**: If you expected `NA` values to be silently ignored, use `na.rm = TRUE` or `na_strategy = "remove"`. |
| **System Dependency (TBB)** | Required linking to and headers for Intel TBB through `RcppParallel`. | **Removed**: Replaced by `RcppThread` (C++17 STL). TBB-specific environment flags no longer have any effect. | ℹ️ **Positive**: No breaking-change risk; TBB environment variables no longer need to be configured. |

---

## 3. `tsmp`: What Changed (0.4.16 $\to$ 0.5.0)

### 3.1 Major Breaking Changes & Critical Removals

> [!WARNING]
> **Pay attention to the dependency and parallelism changes below when updating legacy scripts.**

#### 1. Complete Removal of `doSNOW` and `foreach`
* **Before**: `mstomp_par()` required explicitly creating and registering socket clusters with `parallel::makeCluster()` and `doSNOW::registerDoSNOW()`.
* **Now**: The `doSNOW` and `foreach` packages have been **completely removed** from `tsmp` dependencies. Parallelism is now handled by native in-memory multithreading through `n_workers`.
* **Impact**: Scripts that passed cluster objects or called `registerDoSNOW()` before `mstomp_par()` no longer need to do so. `mstomp_par()` now accepts `n_workers = 2L` directly (or any desired number of threads).

#### 2. Removal of the `audio` Package
* **Before**: Functions played beeps and sound alerts when lengthy computations finished, using the `audio` package.
* **Now**: The `audio` package and all sound calls have been **completely removed**. Scripts that rely on sound side effects should use operating system or R notifications instead, if desired.

#### 3. Invalid Profile Indices Standardized to `NA_integer_`
* **Before**: In `MatrixProfile` or `MultiMatrixProfile` objects, indices with no neighbor could appear as `-1`.
* **Now**: They are standardized to `NA_integer_`.

#### 4. Package Is Now Pure R (`NeedsCompilation: no`)
* **Impact**: `tsmp` now requires `Imports: matrixprofiler (>= 0.1.12)`. C++ build tools are no longer required to install `tsmp` from source.

---

### 3.2 Signature Changes and New Capabilities

#### 1. `find_snippet()`
* **New signature**:
  `find_snippet(data, s_size, n_snippets = 2L, window_size = s_size / 2L, n_workers = 1L)`
* **What's new**: Added the `n_workers` parameter for parallel MPdist computation.
* **Bug fix**: Fixed the longstanding `data is too long` error caused by `matrix(0, nrow = 0)` when the time series length was exactly divisible by `s_size`.

#### 2. `tsmp()` and `compute()`
* **New mode**: `tsmp(..., mode = "mpx")` is now a first-class algorithm, eliminating the need for FFT.
* **Missing-value handling**: Added the `na_strategy = c("monolithic", "segmented", "segmented_native")` argument to `tsmp()`.
* **Unification**: `compute()` now maps transparently to MPF structures, producing standard S3 `MatrixProfile` and `PMP` objects.

#### 3. `fluss_cac()`
* **Optimization**: Internally vectorized using `tabulate()`. The signature and return value remain 100% unchanged, while practical complexity in R has been drastically reduced by eliminating slow loops.

#### 4. Concurrency Safeguards
* Validation of `n_jobs` and `n_workers` in `analyze()`, `pmp()`, and `find_snippet()` now uses strict `checkmate::qassert("X1[1,)")` rules, preventing crashes when `parallel::detectCores()` returns `NA` in Docker containers or restricted environments.

---

## 4. Migration Guide: Practical Code Examples

### Case 1: Multithreading in Multivariate STOMP (`mstomp`)

#### ❌ How it worked in `tsmp` 0.4.x:
```r
library(tsmp)
library(doSNOW)
library(parallel)

cl <- makeCluster(4)
registerDoSNOW(cl)

# Required exporting variables to the socket nodes
mp <- mstomp_par(multivariate_data, window_size = 50)

stopCluster(cl)
```

#### ✅ How to do it in `tsmp` 0.5.0:
```r
library(tsmp)

# Direct native parallelism using shared C++ threads (no socket overhead!)
mp <- mstomp(multivariate_data, window_size = 50, n_workers = 4)

# Or use the equivalent signature:
mp <- mstomp_par(multivariate_data, window_size = 50, n_workers = 4)
```

---

### Case 2: Checking for Invalid Matrix Profile Indices

#### ❌ Legacy code:
```r
# Finding indices that had no nearest neighbor
invalid_indices <- which(mp$pi == -1 | mp$pi <= 0)
```

#### ✅ In the new version (in both `tsmp` and `matrixprofiler`):
```r
# Now uses R's idiomatic missing-data convention
invalid_indices <- which(is.na(mp$pi))
```

---

### Case 3: Handling NAs in `matrixprofiler` Moving Window Functions

#### ❌ In `matrixprofiler` 0.1.10:
```r
# If the series contained NAs, mov_mean could produce cascading NaNs or failures
moving_average <- mov_mean(series_with_nas, window_size = 20)
```

#### ✅ How to control this explicitly in `matrixprofiler` 0.1.12:
```r
# 1. Default: Windows that overlap an NA return NA_real_; valid windows remain precise
moving_average <- mov_mean(series_with_nas, window_size = 20, na_strategy = "window_na")

# 2. Ignore NAs when calculating each window
moving_average_without_na <- mov_mean(series_with_nas, window_size = 20, na.rm = TRUE)

# 3. Raise an error immediately if NAs are present
# moving_average_error <- mov_mean(series_with_nas, window_size = 20, na_strategy = "error")
```

---

### Case 4: Time Series Snippets with Parallelism

#### ❌ In `tsmp` 0.4.x:
```r
# Sequential execution, and prone to errors if length(ts) % s_size == 0
snips <- find_snippet(series, s_size = 100, n_snippets = 3)
```

#### ✅ In `tsmp` 0.5.0:
```r
# Multithreading is passed directly to the MPdist C++ engine
snips <- find_snippet(series, s_size = 100, n_snippets = 3, n_workers = 4)
```

---

## 5. Executive Summary for Release Notes / User Announcement

| Feature | Status in 0.4.x / 0.1.10 | Status in 0.5.0 / 0.1.12 | User Action Required |
| :--- | :--- | :--- | :--- |
| **`doSNOW`/`foreach` clusters** | Required for parallel computation | **Removed** | Remove `makeCluster()` and use `n_workers = N`. |
| **Audio alerts (`audio`)** | Present | **Removed** | No technical action; sound alerts have been discontinued. |
| **Invalid profile index** | `-1` / `-1L` | **`NA_integer_`** | Replace `x == -1` with `is.na(x)`. |
| **Moving windows with NA** | Undefined behavior | `na_strategy = "window_na"` | Add `na.rm = TRUE` if you want to ignore NAs. |
| **C++ compiler for `tsmp`** | Required when installing from source | **Not required** (Pure R) | None; installation is now much faster. |
| **MPX in `tsmp()`** | Available only through a wrapper | Native mode `mode = "mpx"` | Use `tsmp(data, w, mode = "mpx")`. |
| **MPdist & snippets** | Pure R / Slow | Native C++ and multicore | Use the `n_workers = N` parameter. |
