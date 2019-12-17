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
analyze <- function(ts, query = None, windows = NULL, sample_pct = 1.0, threshold = 0.98, n_jobs = -1) {

}
