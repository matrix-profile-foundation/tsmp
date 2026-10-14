#' Multivariate STOMP algorithm
#'
#' Computes the Matrix Profile and Profile Index for Multivariate Time Series.
#'
#' @details
#' The Matrix Profile, has the potential to revolutionize time series data mining because of its
#' generality, versatility, simplicity and scalability. In particular it has implications for time
#' series motif discovery, time series joins, shapelet discovery (classification), density
#' estimation, semantic segmentation, visualization, rule discovery, clustering etc. The MSTOMP
#' computes the Matrix Profile and Profile Index for Multivariate Time Series that is meaningful for
#' multidimensional MOTIF discovery. It uses the STOMP algorithm that is faster than STAMP but lacks
#' its anytime property.
#'
#' Although this functions handles Multivariate Time Series, it can also be used to handle
#' Univariate Time Series. `verbose` changes how much information is printed by this function; `0`
#' means nothing, `1` means text, `2` adds the progress bar.
#'
#' @param data a `matrix` of `numeric`, where each column is a time series. Accepts `vector` (see
#'   details), `list` and `data.frame` too.
#' @param window_size an `int` with the size of the sliding window.
#' @param exclusion_zone a `numeric`. Size of the exclusion zone, based on window size (default is
#'   `1/2`).
#' @param verbose an `int`. See details. (Default is `2`).
#' @param must_dim an `int` or `vector` of which dimensions to forcibly include (default is `NULL`).
#' @param exc_dim an `int` or `vector` of which dimensions to exclude (default is `NULL`).
#' @param n_workers an `int`. Number of threads for parallel execution (default is `1`).
#'
#' @return Returns a `MultiMatrixProfile` object, a `list` with the matrix profile `mp`, profile index `pi`
#'   left and right matrix profile `lmp`, `rmp` and profile index `lpi`, `rpi`, window size `w`,
#'   number of dimensions `n_dim`, exclusion zone `ez`, must dimensions `must` and excluded dimensions `exc`.
#'
#'   If the input has only one dimension, returns the same as [stomp()].
#'
#' @export
#'
#' @describeIn mstomp Single thread version.
#'
#' @family matrix profile computations
#'
#' @references * Yeh CM, Kavantzas N, Keogh E. Matrix Profile VI : Meaningful Multidimensional Motif
#'   Discovery.
#' @references * Zhu Y, Imamura M, Nikovski D, Keogh E. Matrix Profile VII: Time Series Chains: A
#'   New Primitive for Time Series Data Mining. Knowl Inf Syst. 2018 Jun 2;1-27.
#' @references Website: <https://sites.google.com/view/mstamp/>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#' # using all dimensions
#' mp <- mstomp(mp_toy_data$data[1:150, ], 30, verbose = 0)
#' \donttest{
#' # using threads
#' mp <- mstomp_par(mp_toy_data$data[1:150, ], 30, verbose = 0)
#'
#' # force using dimensions 1 and 2
#' mp <- mstomp(mp_toy_data$data[1:200, ], 30, must_dim = c(1, 2))
#' # exclude dimensions 2 and 3
#' mp2 <- mstomp(mp_toy_data$data[1:200, ], 30, exc_dim = c(2, 3))
#' }
#'
mstomp <- function(data, window_size, exclusion_zone = getOption("tsmp.exclusion_zone", 1 / 2),
                   verbose = getOption("tsmp.verbose", 2),
                   must_dim = NULL, exc_dim = NULL, n_workers = 1L) {
  tictac <- Sys.time()
  progress <- isTRUE(verbose >= 2)

  res <- matrixprofiler::mstomp(
    data = data,
    window_size = window_size,
    exclusion_zone = exclusion_zone,
    must_dim = must_dim,
    exc_dim = exc_dim,
    n_workers = n_workers,
    progress = progress
  )

  tictac <- Sys.time() - tictac
  if (verbose > 0) {
    message(sprintf("Finished in %.2f %s", tictac, units(tictac)))
  }

  if (res$n_dim > 1L) {
    obj <- list(
      mp = res$matrix_profile,
      pi = res$profile_index,
      rmp = res$right_matrix_profile,
      rpi = res$right_profile_index,
      lmp = res$left_matrix_profile,
      lpi = res$left_profile_index,
      w = res$w,
      ez = res$ez,
      n_dim = res$n_dim,
      must = res$must,
      exc = res$exc
    )
    class(obj) <- "MultiMatrixProfile"
    attr(obj, "join") <- FALSE
  } else {
    obj <- list(
      mp = res$matrix_profile,
      pi = res$profile_index,
      rmp = res$right_matrix_profile,
      rpi = res$right_profile_index,
      lmp = res$left_matrix_profile,
      lpi = res$left_profile_index,
      w = res$w,
      ez = res$ez
    )
    class(obj) <- "MatrixProfile"
    attr(obj, "join") <- FALSE
  }

  return(obj)
}
