#' Contrast Profile
#'
#' Computes the contrast profile of two (classes of) time series.
#'
#' @param negative_data Required. Any 1-dimension series of numbers (`matrix`, `vector`, `ts` etc.) where the pattern is not present
#' @param positive_data Required. Any 1-dimension series of numbers (`matrix`, `vector`, `ts` etc.) where the pattern is present
#' @param window_size Required. An integer defining the rolling window size.
#' @param positive_matrix Optional. A precomputed self-similar matrix profile of the positive data.
#' @param exclusion_zone A numeric. Defines the size of the area around the rolling window that will be ignored to avoid
#'   trivial matches. Default is `0.5`, i.e., half of the `window_size`.
#' @param distance A string. Currently accepts `euclidean` and `pearson`. Defaults to `euclidean`.
#' @param n_workers An integer. The number of threads using for computing. Defaults to `1`.
#' @param progress A logical. If `TRUE` (the default) will show a progress bar. Useful for long computations.
#'
#' @details ## Constrast Profile
#' This algorithm returns the contrast profile of two time series, which shows the position of patters that are similar in the
#' positive data, but at the same time very dissimilar in the negative data.  In other words, this means that such a pattern
#' represents well positive data and may be taken as a "signature" of that class. More information can be found in the references.
#'
#' @return Returns a `list` with the `contrast_profile`, `plato`, `plato_nn`, `plato_idx`, `plato_nn_idx`, `w`, `ez`, `euclidean`
#' values
#'
#' @references * R. Mercer, S. Alaee, A. Abdoli, S. Singh, A. Murillo and E. Keogh, "Matrix Profile XXIII: Contrast Profile:
#'   A Novel Time Series Primitive that Allows Real World Classification," 2021 IEEE International Conference on Data Mining (ICDM), 2021,
#'   pp. 1240-1245, doi: 10.1109/ICDM51629.2021.00151.
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @export
#' @examples
#' cp <- contrast(motifs_discords_small, rev(motifs_discords_small), 50)
#'
contrast <- function(negative_data, positive_data, window_size, positive_matrix = NULL,
                     exclusion_zone = getOption("tsmp.exclusion_zone", 1 / 2), distance = c("euclidean", "pearson"),
                     n_workers = 1L, progress = TRUE, verbose = getOption("tsmp.verbose", 2)) {
  "!!!DEBUG Parsing Arguments"
  negative_data <- as.numeric(negative_data)
  checkmate::qassert(negative_data, "N+")
  positive_data <- as.numeric(positive_data)
  checkmate::qassert(positive_data, "N+")
  window_size <- as.integer(checkmate::qassert(
    window_size,
    "X+"
  ))
  checkmate::qassert(positive_matrix, c("0", "L"))

  if (!is.null(positive_matrix)) {
    if (!("MatrixProfile" %in% class(positive_matrix))) {
      cli::cli_abort("`positive_matrix` argument must be an object of class `MatrixProfile`.")
    }
  } else {
    positive_matrix <- list()
  }

  checkmate::qassert(exclusion_zone, "N+")
  distance <- match.arg(distance)
  if (distance == "euclidean") {
    dist <- TRUE
  } else {
    dist <- FALSE
  }
  n_workers <- as.integer(checkmate::qassert(n_workers, "X+"))
  checkmate::qassert(progress, "B+")

  ez <- exclusion_zone
  result <- NULL

  smaller_size <- min(length(negative_data), length(positive_data))
  if (window_size > ceiling(smaller_size / 2)) {
    cli::cli_abort("The smaller time series is too short relative to desired window size.")
  }

  "!DEBUG Register anytime exit point"
  on.exit(if (is.null(result)) {
    return(invisible(NULL))
  } else {
    return(result)
  }, TRUE)

  "!DEBUG Computation"
  tryCatch(
    {
      "!DEBUG n_workers = `n_workers`"
      if (n_workers > 1) {
        p <- RcppParallel::defaultNumThreads()
        n_workers <- min(n_workers, p)
        RcppParallel::setThreadOptions(numThreads = n_workers)
        result <- matrixprofiler:::contrast_profile_rcpp(negative_data, positive_data,
          window_size, positive_matrix,
          ez = ez, s_size = 1,
          n_workers = as.integer(n_workers), as.logical(dist),
          as.logical(progress)
        )
        RcppParallel::setThreadOptions(numThreads = p)
      } else {
        result <- matrixprofiler:::contrast_profile_rcpp(negative_data, positive_data,
          window_size, positive_matrix,
          ez = ez, s_size = 1,
          n_workers = 1L, as.logical(dist), as.logical(progress)
        )
      }
    },
    error = print
  )
}
