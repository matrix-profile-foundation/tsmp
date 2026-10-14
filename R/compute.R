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
  data_ref <- convert_data(ts)
  query_ref <- convert_data(query)

  # Parse arguments ---------------------------------
  checkmate::qassert(data_ref, "N+")
  windows <- checkmate::qassert(windows, c("0", "X+[4,)"))
  checkmate::qassert(query_ref, c("0", "N>=4"))
  checkmate::qassert(sample_pct, "N1(0,1]")
  checkmate::qassert(threshold, c("0", "N1(0,1]"))
  n_jobs <- as.integer(checkmate::qassert(n_jobs, "X1[1,)"))

  if (length(windows) == 1) {
    ## Matrix Profile (delegated to tsmp) ======
    if (is.null(query)) {
      result <- tsmp(ts,
        window_size = windows,
        mode = "mpx",
        s_size = sample_pct,
        n_workers = n_jobs,
        verbose = 0,
        .keep_data = FALSE
      )
    } else {
      result <- tsmp(ts, query,
        window_size = windows,
        mode = "mpx",
        s_size = sample_pct,
        n_workers = n_jobs,
        verbose = 0,
        .keep_data = FALSE
      )
    }

    result$sample_pct <- sample_pct
    result$data <- list(
      ts = ts,
      query = query
    )
    attr(result, "metric") <- "euclidean"
    attr(result, "algorithm") <- "mpx"
    attr(result, "join") <- !is.null(query)
    class(result) <- "MatrixProfile"
    return(invisible(result))
  } else {
    ## Pan Matrix Profile =========================
    if (!is.null(threshold)) {
      # when a threshold is passed, we compute the upper bound
      res <- pmp_upper_bound(data = data_ref, threshold = threshold, n_workers = n_jobs)
    }

    if (is.null(windows)) {
      # when no windows are passed, create an array from 10 to upper bound or half ts size
      windows <- seq.int(from = 10, to = min(length(data_ref) / 2, res$upper_window), length.out = 20)
    } else {
      # otherwise, remove windows that are above upper bound or half ts size
      windows <- windows[windows <= min(length(data_ref) / 2, res$upper_window)]
    }

    windows <- floor(windows)
    result <- pmp(data = data_ref, window_sizes = windows, plot = FALSE, pmp_obj = res, n_workers = n_jobs, verbose = 0)
    result$sample_pct <- sample_pct
    result$data <- list(
      ts = ts,
      query = query
    )
    attr(result, "join") <- FALSE
    attr(result, "metric") <- "euclidean"
    attr(result, "algorithm") <- "pmp"
    class(result) <- "PMP"
    return(invisible(result))
  }
}
