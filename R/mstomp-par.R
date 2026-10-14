#' Multivariate STOMP algorithm Parallel version
#'
#' @param n_workers an `int`. Number of workers for parallel. (Default is `2`).
#'
#' @export
#'
#' @describeIn mstomp Parallel version.
#'
mstomp_par <- function(data, window_size, exclusion_zone = getOption("tsmp.exclusion_zone", 1 / 2),
                       verbose = getOption("tsmp.verbose", 2),
                       must_dim = NULL, exc_dim = NULL, n_workers = 2L) {
  mstomp(
    data = data,
    window_size = window_size,
    exclusion_zone = exclusion_zone,
    verbose = verbose,
    must_dim = must_dim,
    exc_dim = exc_dim,
    n_workers = n_workers
  )
}
