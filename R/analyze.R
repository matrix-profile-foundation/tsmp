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
#' 1. Compute the upper bound when a threshold is provided
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
analyze <- function(ts, windows = NULL, query = NULL, sample_pct = 1.0, threshold=0.98, n_jobs = 1L) {

  # Parse arguments ---------------------------------
  checkmate::qassert(ts, "N+")
  windows <- checkmate::qassert(windows, c("0", "X+"))
  checkmate::qassert(query, c("0", "N>=4"))
  checkmate::qassert(sample_pct, "N1(0,1]")
  checkmate::qassert(threshold, c("0", "N1(0,1]"))
  n_jobs <- as.integer(checkmate::qassert(n_jobs, paste0("X1[1,", parallel::detectCores(), "]")))


  # Compute --------------------------------------
  result <- compute(ts, windows, query, sample_pct, threshold, n_jobs)

  # Explore --------------------------------------
  # extract top motifs
  result <- find_motif(result, n_motifs = 5) %T>% plot()

  # extract top discords
  result <- find_discord(result, n_discords = 5) %T>% plot()

  # Visualize --------------------------------------
  figures <- NULL
  #figures <- visualize(result)

  return(list(profile = result, figures = figures))
}
