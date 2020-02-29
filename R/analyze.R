#' Runs an appropriate workflow based on the parameters passed in.
#'
#' The goal of this function is to compute all fundamental algorithms on the provided
#' time series data. See details for more information.
#'
#' @inheritParams compute
#'
#' @details
#' For now the following is computed:
#'
#' 1. Matrix Profile - exact or approximate based on `sample_pct` given that a single `windows` is provided. By default
#' is the exact algorithm;
#' 2. Top 3 Motifs;
#' 3. Top 3 Discords;
#' 4. Plot Matrix Profile, Motifs and Discords.
#'
#' When `windows` is not provided or more than a single window is provided,
#' the Pan-Matrix Profile is computed:
#'
#' 1. Compute the upper bound when a `threshold` is provided (it is, by default);
#' 2. Compute Pan-Matrix Profile for all `windows` provided, below the upper bound, or a default range when no `windows`
#' is provided;
#' 3. Top Motifs;
#' 4. Top Discords;
#' 5. Plot Pan-Matrix Profile, motifs and discords.
#'
#' @return
#' The appropriate Matrix Profile or Pan-Matrix Profile profile object and also plots the graphics.
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
#' result <- analyze(mp_toy_data$data[, 1], 80)
#' \dontrun{
#' # Pan Matrix Profile
#' result <- analyze(mp_toy_data$data[, 1])
#' }
analyze <- function(ts, windows = NULL, query = NULL, sample_pct = 1.0, threshold = 0.98, n_jobs = 1L) {

  # Parse arguments ---------------------------------
  checkmate::qassert(ts, "N+")
  windows <- checkmate::qassert(windows, c("0", "X+[4,)"))
  checkmate::qassert(query, c("0", "N>=4"))
  checkmate::qassert(sample_pct, "N1(0,1]")
  checkmate::qassert(threshold, c("0", "N1(0,1]"))
  n_jobs <- as.integer(checkmate::qassert(n_jobs, paste0("X1[1,", parallel::detectCores(), "]")))

  result <- NULL

  # Register the anytime exit point
  on.exit(
    {
      return(invisible(result))
    },
    TRUE
  )

  # Compute --------------------------------------
  result <- compute(ts, windows, query, sample_pct, threshold, n_jobs) %T>% visualize()

  # Explore --------------------------------------
  # extract top motifs
  result <- motifs(result, k = 3L) %T>% visualize()

  # extract top discords
  result <- discords(result, k = 3L) %T>% visualize()
}
