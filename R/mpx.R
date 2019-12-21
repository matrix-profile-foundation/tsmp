
#' Fast implementation of MP and MPI for internal purposes, without FFT
#'
#' @param data data
#' @param window_size window size
#' @param query query
#' @param idx calc and return indexes?
#' @param dist
#' @param n_workers
#'
#' @return Returns MP and MPI
#' @keywords internal
#' @noRd

mpx <- function(data, window_size, query = NULL, idx = TRUE, dist = c("euclidean", "pearson"), n_workers = 1) {

  # Parse arguments ---------------------------------
  minlag <- floor(window_size / 4)
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

  result <- NULL

  # Register anytime exit point
  on.exit(return({
    result
  }), TRUE)

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
    tryCatch(
      {
        if (n_workers > 1) {
          # TODO: AB multithread
          message("AB Multi-thread not yet implemented, using 1 thread instead.")
          p <- RcppParallel::defaultNumThreads()
          n_workers <- min(n_workers, p)
          RcppParallel::setThreadOptions(numThreads = n_workers)
          result <- mpxab_rcpp(
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
