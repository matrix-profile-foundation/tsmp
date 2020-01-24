#' Search for Motifs
#'
#' @param profile the PMP or MP data structure
#' @param exclusion_zone an `int`. Number of values to exclude on both sides of the motif to avoid trivial matches.
#' Defaults to the exclusion zone used to compute the MP/PMP which is found in the profile data structure.
#' @param k an `int`. Number of motifs to find. (Default is `3`).
#' @param neighbor_count an `int`. Number of neighbors to find. (Default is `3`).
#' @param radius an `int`. Set a threshold to exclude matching neighbors with distance > current
#' discord distance * `radius`. (Default is `3`).
#'
#' @name motifs
#' @export

motifs <- function(profile, exclusion_zone = profile$ez, k = 3L, neighbor_count = 10L, radius = 3) {
  find_motif(profile, n_motifs = k, n_neighbors = neighbor_count, radius = radius, exclusion_zone = exclusion_zone)
}
