
#' Fast implementation of MP and MPI for internal purposes, without FFT
#'
#' @param data a `matrix` or a `vector`. The time series to analyze.
#' @param window_size window size
#' @param query query
#' @param idx  compute the profile indexes?
#' @param dist distance measure, Euclidean or Pearson?
#' @param n_workers threads for multi-threading
#'
#' @return Returns MP and MPI
#' @export

mpx <- function(data, window_size, query = NULL, idx = TRUE, dist = c("euclidean", "pearson"), n_workers = 1) {

  # Parse arguments ---------------------------------
  minlag <- floor(window_size / 2)
  dist <- match.arg(dist)
  checkmate::qassert(data, "N+")
  window_size <- as.integer(checkmate::qassert(window_size, "X+"))
  n_workers <- as.integer(checkmate::qassert(n_workers, "X+"))
  checkmate::qassert(query, c("0", "N>=4"))

  if (dist == "euclidean") {
    dist <- TRUE
  } else {
    dist <- FALSE
  }

  ez <- getOption("tsmp.exclusion_zone", 1 / 2) # minlag is the exclusion zone
  result <- NULL

  # Register anytime exit point
  on.exit(
    if (is.null(result)) {
      return(invisible(NULL))
    } else {
      result$ez <- ez
      return(result)
    },
    TRUE
  )

  # Computation ------------------------------------
  if (is.null(query)) {
    ## Self-Join ====================================
    tryCatch(
      {
        if (n_workers > 1) {
          p <- RcppParallel::defaultNumThreads()
          n_workers <- min(n_workers, p)
          RcppParallel::setThreadOptions(numThreads = n_workers)
          result <- mpx_rcpp_parallel(
            data,
            window_size,
            as.integer(minlag),
            as.logical(idx),
            as.logical(dist)
          )
          RcppParallel::setThreadOptions(numThreads = p)
        } else {
          result <- mpx_rcpp(
            data,
            window_size,
            as.integer(minlag),
            as.logical(idx),
            as.logical(dist)
          )
        }
      },
      error = print
    )
  } else {
    ## AB-Join ====================================
    ez <- 0

    tryCatch(
      {
        if (n_workers > 1) {
          p <- RcppParallel::defaultNumThreads()
          n_workers <- min(n_workers, p)
          RcppParallel::setThreadOptions(numThreads = n_workers)
          result <- mpxab_rcpp_parallel(
            data,
            query,
            window_size,
            as.logical(idx),
            as.logical(dist)
          )
          RcppParallel::setThreadOptions(numThreads = p)
        } else {
          result <- mpxab_rcpp(
            data,
            query,
            window_size,
            as.logical(idx),
            as.logical(dist)
          )
        }
      },
      error = print
    )
  }
}
