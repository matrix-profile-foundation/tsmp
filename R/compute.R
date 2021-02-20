#' Computes the Matrix Profile or Pan-Matrix Profile
#'
#' Main API Function
#'
#' Computes the exact or approximate Matrix Profile based on the sample percent
#' specified. Currently, MPX and SCRIMP++ are used for the exact and
#' approximate algorithms respectively. See details for more information about the arguments
#' combinations.
#'
#' @param ts a `matrix` or a `vector`. The time series to analyze.
#' @param query a `matrix` or a `vector`. Optional The query to analyze. Note that when computing the Pan-Matrix Profile
#'  the query is ignored!
#' @param windows an `int` or a `vector`. The window(s) to compute the Matrix Profile. Note that it may be an `int`
#' for a single matrix profile computation or a `vector` of `int` for computing the Pan-Matrix Profile.
#' @param sample_pct a `numeric`. A number between 0 and 1 representing how many samples to compute for
#' the Matrix Profile or Pan-Matrix Profile. When it is 1, the exact algorithm is used. (default is `1.0`).
#' @param threshold a `numeric`. Correlation threshold. See details.  (Default is `0.98`).
#' @param n_jobs an `int`. The number of cpu cores to use when computing the MatrixProfile. (default is `1`).
#'
#' @details
#'
#' When a single `windows` is given, the Matrix Profile is computed. If a `query` is provided, AB join is computed.
#' Otherwise the self-join is computed.
#' When multiple `windows` or none are given, the Pan-Matrix Profile is computed. If a `threshold` is set (it is,
#' by default), the upper bound will be computed and the given `windows` or a default range (when no `windows`), below
#' the upper bound will be computed.
#'
#' @return
#' The profile computed.
#'
#' @export
#'
#' @family Main API
#'
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#'
#' # Matrix Profile
#' result <- compute(mp_toy_data$data[, 1], 80)
#' \donttest{
#' # Pan-Matrix Profile
#' result <- compute(mp_toy_data$data[, 1])
#' }
compute <- function(ts, windows = NULL, query = NULL, sample_pct = 1.0, threshold = 0.98, n_jobs = 1L) {

  # Check time series classes -----------------------
  ## ts ----
  ts_class <- class(ts)
  if (ts_class %in% c("data.frame", "list")) {
    stop(glue::glue("Time Series cannot be `{ts_class}` and must have a single dimension."))
  }
  ts_attr <- attributes(ts)
  ts_names <- names(ts)
  ts <- as.vector(ts)

  if (!is.null(query)) {
    query_class <- class(query)
    if (query_class %in% c("data.frame", "list")) {
      stop(glue::glue("Query cannot be `{query_class}` and must have a single dimension."))
    }
    query_attr <- attributes(query)
    query_names <- names(query)
    query <- as.vector(query)
  }

  # Parse arguments ---------------------------------
  checkmate::qassert(ts, "N+")
  windows <- checkmate::qassert(windows, c("0", "X+[4,)"))
  checkmate::qassert(query, c("0", "N>=4"))
  checkmate::qassert(sample_pct, "N1(0,1]")
  checkmate::qassert(threshold, c("0", "N1(0,1]"))
  n_jobs <- as.integer(checkmate::qassert(n_jobs, paste0("X1[1,", 4, "]")))

  res <- NULL
  algorithm <- NULL
  join <- FALSE
  metric <- "euclidean"

  on.exit({
    if (!is.null(res)) {
      # Build compute object --------------------------

      # Main fields, easily accessible by the user
      if (length(windows) == 1) {
        result <- list(
          mp = res$matrix_profile,
          pi = res$profile_index,
          # mpb = NULL,
          # pib = NULL,
          rmp = NULL,
          rpi = NULL,
          lmp = NULL,
          lpi = NULL,
          w = windows,
          ez = res$ez
        )
        class(result) <- "MatrixProfile"
      } else {
        result <- res
        class(result) <- "PMP"
      }
      result$partial <- res$partial
      result$sample_pct <- sample_pct

      class(ts) <- ts_class
      attributes(ts) <- ts_attr
      names(ts) <- ts_names

      if (!is.null(query)) {
        class(query) <- query_class
        attributes(query) <- query_attr
        names(query) <- query_names
      }

      result$data <- list(
        ts = ts,
        query = query
      )

      # Attributes
      attr(result, "join") <- join
      attr(result, "metric") <- metric
      attr(result, "algorithm") <- algorithm

      # End ---------------------------------
      return(invisible(result))
    }
  })

  # Start ---------------------------------
  if (length(windows) == 1) {
    ## Matrix Profile =========================

    if (is.null(query)) {
      ### Self-join #############################

      if (sample_pct >= 1) {
        algorithm <- "mpx"
        res <- matrixprofiler::mpx(data = ts, window_size = windows, n_workers = n_jobs, progress = TRUE)
      } else {
        algorithm <- "mpx" # former was scrimp
        res <- matrixprofiler::mpx(
          data = ts, window_size = windows, s_size = sample_pct,
          n_workers = n_jobs, progress = TRUE
        )
      }
    } else {
      ### AB join #############################
      join <- TRUE
      if (sample_pct >= 1) {
        algorithm <- "mpx"
        res <- matrixprofiler::mpx(data = ts, query = query, window_size = windows, n_workers = n_jobs, progress = TRUE)
      } else {
        algorithm <- "scrimp"
        res <- matrixprofiler::mpx(
          data = ts, window_size = windows, s_size = sample_pct,
          n_workers = n_jobs, progress = TRUE
        )
      }
    }
  } else {
    ## Pan Matrix Profile =========================
    algorithm <- "pmp"

    if (!is.null(threshold)) {
      # when a threshold is passed, we compute the upper bound
      res <- pmp_upper_bound(data = ts, threshold = threshold, n_workers = n_jobs)
    }

    if (is.null(windows)) {
      # when no windows are passed, create an array from 10 to upper bound or half ts size
      windows <- seq.int(from = 10, to = min(length(ts) / 2, res$upper_window), length.out = 20)
    } else {
      # otherwise, remove windows that are above upper bound or half ts size
      windows <- windows[windows <= min(length(ts) / 2, res$upper_window)]
    }

    windows <- floor(windows)
    res <- pmp(data = ts, window_sizes = windows, plot = FALSE, pmp_obj = res, n_workers = n_jobs)
  }
}
