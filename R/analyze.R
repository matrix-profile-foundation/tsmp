#' Runs an appropriate workflow based on the parameters passed in.
#'
#' The goal of this function is to compute all fundamental algorithms on the provided
#' time series data. For now the following is computed:
#'
#' 1. Matrix Profile - exact or approximate based on sample_pct given that a
#' window is provided. By default is the exact algorithm.
#' 2. Top Motifs - The top 3 motifs are found.
#' 3. Top Discords - The top 3 discords are found.
#' 4. Plot MP, Motifs and Discords
#'
#' When a window is not provided or more than a single window is provided,
#' the PMP is computed:
#'
#' 1. Compute UPPER window when no window(s) is provided
#' 2. Compute PMP for all windows
#' 3. Top Motifs
#' 4. Top Discords
#' 5. Plot PMP, motifs and discords.
#'
#' @param ts
#' @param query
#' @param windows
#' @param sample_pct
#' @param threshold
#' @param n_jobs
#'
#' @return
#' tuple : (profile, figures) The appropriate PMP or MP profile object and associated figures.
#' @export
#'
#' @examples
analyze <- function(ts, query = NULL, windows = NULL, sample_pct = 1.0, threshold = 0.98, n_jobs = 1) {

  # Parse arguments ---------------------------------
  checkmate::qassert(ts, "N+")
  windows <- checkmate::qassert(windows, c("0", "X+"))
  checkmate::qassert(query, c("0", "N>=4"))
  checkmate::qassert(sample_pct, "N1(0,1]")
  checkmate::qassert(threshold, "N1(0,1]")
  n_jobs <- as.integer(checkmate::qassert(n_jobs, paste0("X1[1,", parallel::detectCores(), "]")))


  # Compute --------------------------------------
  # TODO: I think all this could be inside compute()
  if (is.null(windows) || length(windows) > 1) {
    # PMP when no window provided or several windows
    # when a threshold is passed, we compute the upper window

    result <- NULL

    if (is.null(windows)) {
      result <- tsmp::pmp_upper_bound(data = ts, threshold = threshold, n_workers = n_jobs)
      windows <- seq.int(from = 10, to = length(ts) / 2, length.out = 20)
    }

    result <- tsmp::pmp(data = ts, window_sizes = windows, plot = FALSE, pmp_obj = result, n_workers = n_jobs)
  } else if (sample_pct < 1) {
    # MP approximate
    result <- tsmp::mpx(data = ts, query = query, window_size = windows, idx = TRUE, dist = metric, n_workers = n_jobs)
  } else {
    # MP exact
    result <- scrimp(ts, window_size = windows, s_size = floor(sample_pct * length(ts))) # n_jobs
  }

  # Explore --------------------------------------
  # extract top motifs
  result <- top_k_motifs(result)

  # extract top discords
  result <- top_k_discords(result)

  # Visualize --------------------------------------
  figures <- visualize(result)

  return(list(profile = result, figures = figures))
}
