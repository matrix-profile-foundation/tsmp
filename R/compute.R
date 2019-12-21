#' Computes the exact or approximate MatrixProfile
#'
#' Computes the exact or approximate MatrixProfile based on the sample percent
#' specified. Currently, MPX and SCRIMP++ is used for the exact and
#' approximate algorithms respectively. When multiple windows are passed, the
#' Pan-MatrixProfile is computed and returned.
#'
#' @param ts a `matrix` or a `vector`. The time series to analyze.
#' @param windows an `int` or a `vector`. The window(s) to compute the MatrixProfile. Note that it may be an `int`
#' for a single matrix profile computation or a `vector` of `int` for computing the pan matrix profile (PMP).
#' @param query a `matrix` or a `vector`. Optional The query to analyze. Note that when computing the PMP the query
#' is ignored!
#' @param sample_pct a `numeric`. A number between 0 and 1 representing how many samples to compute for
#' the MP or PMP. When it is 1, the exact algorithm is used. (default is `1`).
#' @param n_jobs an `int`. The number of cpu cores to use when computing the MP. (default is `1`).
#'
#' @return
#' The profile computed.

#' @export
#'
#' @examples
compute <- function(ts, windows = NULL, query = NULL, sample_pct = 1, threshold=0.98, n_jobs = 1) {

  # Parse arguments ---------------------------------
  checkmate::qassert(ts, "N+")
  windows <- as.integer(checkmate::qassert(windows, "X+"))
  checkmate::qassert(query, c("0", "N>=4"))
  checkmate::qassert(sample_pct, "N1(0,1]")
  n_jobs <- as.integer(checkmate::qassert(n_jobs, paste0("X1[1,", parallel::detectCores(), "]")))

  res <- NULL
  algorithm <- NULL
  join <- FALSE
  metric <- "euclidean"

  # Start ---------------------------------
  if (length(windows) == 1) {
    ## Matrix Profile =========================

    if (is.null(query)) {
      ### Self-join #############################

      if (sample_pct >= 1) {
        res <- tsmp::mpx(data = ts, window_size = windows, idx = TRUE, dist = metric, n_workers = n_jobs)
        algorithm <- "mpx"
      } else {
        res <- scrimp(ts, window_size = windows, s_size = floor(sample_pct * length(ts))) # n_jobs
        algorithm <- "scrimp"
      }

    } else {
      ### AB join #############################
      join <- TRUE
      res <- tsmp::mpx(data = ts, query = query, window_size = windows, idx = TRUE, dist = metric, n_workers = n_jobs)
      algorithm <- "mpx"
      # TODO: add scrimp AB-join
    }
  } else {
    ## Pan Matrix Profile =========================

    res <- pmp(ts, window_sizes = windows, plot = FALSE) # n_jobs
    algorithm <- "pmp"
  }

  # Build compute object --------------------------

  result <- list()

  # Main fields, easily accessible by the user
  if (length(windows) == 1) {
    class(result) <- "MatrixProfile"
    result$mp <- res$mp
    result$pi <- res$pi
    result$rmp <- list(NULL)
    result$rpi <- list(NULL)
    result$lmp <- list(NULL)
    result$lpi <- list(NULL)
  } else {
    class(result) <- "PanMatrixProfile"
    result$pmp <- res # TODO
  }
  result$w <- windows
  result$ez <- getOption("tsmp.exclusion_zone", 1 / 2)
  result$data <- list(ts = ts,
                      query = query)

  # Attributes
  attr(result, "join") <- join
  attr(result, "metric") <- metric
  attr(result, "algorithm") <- algorithm

  # End ---------------------------------

  return(result)
}
