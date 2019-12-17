
#' Fast implementation of MP and MPI for internal purposes, without FFT
#'
#' @param data data
#' @param window_size window size
#' @param query
#' @param idx calc and return indexes?
#' @param dist
#'
#' @return Returns MP and MPI
#' @keywords internal
#' @noRd

mpx <- function(data, window_size, query = NULL, idx = FALSE, dist = c("euclidean", "pearson")) {

  # Parse arguments ---------------------------------
  minlag <- floor(window_size / 4)
  dist <- match.arg(dist)
  checkmate::qassert(data, "N+")
  window_size <- as.integer(checkmate::qassert(window_size, "X+"))
  checkmate::qassert(query, c("0", "N>=4"))

  if (dist == "euclidean") {
    dist <- TRUE
  } else {
    dist <- FALSE
  }

  # # Register anytime exit point
  # on.exit(return({
  #   result
  # }), TRUE)

  if (is.null(query)) {
    result <- mpx_rcpp(
      data,
      window_size,
      as.integer(minlag),
      as.logical(idx),
      as.logical(dist)
    )
  } else {
    result <- mpxab_rcpp(
      data,
      query,
      window_size,
      as.logical(idx),
      as.logical(dist)
    )
  }
}
