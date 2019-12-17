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
#' dict : profile
#'     A dict of key data points computed.
#' {
#' 'mp': The matrix profile,
#' 'pi': The matrix profile 1NN indices,
#' 'rmp': The right matrix profile,
#' 'rpi': The right matrix profile 1NN indices,
#' 'lmp': The left matrix profile,
#' 'lpi': The left matrix profile 1NN indices,
#' 'metric': The distance metric computed for the mp,
#' 'w': The window size used to compute the matrix profile,
#' 'ez': The exclusion zone used,
#' 'join': Flag indicating if a similarity join was computed,
#' 'data': {
#' 'ts': Time series data,
#'     'query': Query data if supplied
#'   }
#'   'class': "MatrixProfile"
#'   'algorithm': "mpx"
#' }
#' The profile computed.

#' @export
#'
#' @examples
compute <- function(ts, windows, query = NULL, sample_pct = 1, n_jobs = 1) {

  # Parse arguments ---------------------------------
  checkmate::qassert(ts, "N+")
  windows <- as.integer(checkmate::qassert(windows, "X+"))
  checkmate::qassert(query, c("0", "N>=4"))
  checkmate::qassert(sample_pct, "N1(0,1]")
  n_jobs <- as.integer(checkmate::qassert(n_jobs, paste0("X1[1,", parallel::detectCores(), "]")))

  # Start ---------------------------------
  if (length(windows) == 1) {
    ## Matrix Profile =========================

    if (is.null(query)) {
      ### Self-join #############################

      if (sample_pct >= 1) {
        result <- tsmp(ts, window_size = windows, mode = "stomp", n_workers = n_jobs)
      } else {
        result <- tsmp(ts, window_size = windows, mode = "scrimp", s_size = floor(sample_pct * length(ts)))
      }

    } else {
      ### AB join #############################

      result <- tsmp(ts, query, window_size = windows, mode = "stomp", n_workers = n_jobs)
    }
  } else {
    ## Pan Matrix Profile =========================

    result <- pmp(ts, window_sizes = windows)
  }

  # End ---------------------------------

  return(result)
}
