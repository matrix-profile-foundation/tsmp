#' Fast Low-cost Unipotent Semantic Segmentation (FLUSS)
#'
#' FLUSS is a Domain Agnostic Online Semantic Segmentation that uses the assumption that when few
#' arc are crossing a given index point, means that there is a high probability of semantic change.
#'
#' @details
#' `verbose` changes how much information is printed by this function; `0` means nothing, `1` means
#' text, `2` means text and sound.
#'
#' @param data a `matrix` or a `vector`. Input data.
#' @param window_size an `int`. Size of the sliding window.
#' @param num_segments an `int`. Number of segments to extract. Based on domain knowledge.
#' @param exclusion_zone a `numeric`. Size of the exclusion zone, based on window size (default is
#'   `5`).
#' @param gtruth an `int` or `vector` of `int` with the ground truth index of segments. (Default is
#'   `NULL`).
#' @param profile_index a pre-computed profile index. (Default is `NULL`).
#' @param verbose an `int`. See details. (Default is `2`).
#'
#' @return Returns a list with `segments` (location of semantic changes), `mp` (matrix profile if
#'   computed), `pi` (profile index, input of computed), `cac` corrected arc count.
#' @export
#' @family fluss
#' @references * Gharghabi S, Ding Y, Yeh C-CM, Kamgar K, Ulanova L, Keogh E. Matrix Profile VIII:
#'   Domain Agnostic Online Semantic Segmentation at Superhuman Performance Levels. In: 2017 IEEE
#'   International Conference on Data Mining (ICDM). IEEE; 2017. p. 117–26.
#' @references Website: <https://sites.google.com/site/onlinesemanticsegmentation/>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#' @examples
#' # This is a fast toy example and results are useless. For a complete result, run the code inside
#' #'Not run' section below.
#' data <- mp_fluss_data$tilt_abp$data[1:1000]
#' truth <- 400
#' w <- 10
#' segments <- fluss(data, w, 1, gtruth = truth, verbose = 0)
#' \dontrun{
#' data <- mp_fluss_data$walkjogrun$data
#' w <- mp_fluss_data$walkjogrun$window # 80
#' truth <- mp_fluss_data$walkjogrun$gtruth # 3800 6800
#' nseg <- length(mp_fluss_data$walkjogrun$gtruth) # 2
#' segments <- fluss(data, w, nseg, gtruth = truth)
#' }
fluss <- function(.mp, num_segments, exclusion_zone = NULL) {
  fluss_cac(.mp, exclusion_zone) %>% fluss_extract(num_segments = num_segments, exclusion_zone = exclusion_zone)
}

#' FLUSS - Extract Segments
#'
#' Extract candidate points of semantic changes.
#'
#' @param arc_counts a `matrix` with the corrected arc counts from [fluss_cac()].
#' @param num_segments an `int`. Number of segments to extract. Based on domain knowledge.
#' @param window_size an `int`. Size of the sliding window.
#' @param exclusion_zone a `numeric`. Size of the exclusion zone, based on window size (default is
#'   `5`).
#'
#' @return Returns an `int` or a `vector` of `int` with the location of predicted semantic changes.
#'   The number of locations is not greater than `num_segments`.
#' @export
#' @family fluss
#' @references * Gharghabi S, Ding Y, Yeh C-CM, Kamgar K, Ulanova L, Keogh E. Matrix Profile VIII:
#'   Domain Agnostic Online Semantic Segmentation at Superhuman Performance Levels. In: 2017 IEEE
#'   International Conference on Data Mining (ICDM). IEEE; 2017. p. 117–26.
#' @references Website: <https://sites.google.com/site/onlinesemanticsegmentation/>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#' @examples
#' # This is a fast toy example and results are useless. For a complete result, run the code inside
#' #'Not run' section below.
#' data <- mp_fluss_data$tilt_abp$data[1:1000]
#' w <- 210
#' mp <- stomp(data, window_size = w, verbose = 0)
#' cac <- fluss_cac(mp$pi, w)
#' segments <- fluss_extract(cac, 1, w)
#' \dontrun{
#' data <- mp_fluss_data$walkjogrun$data
#' w <- mp_fluss_data$walkjogrun$window # 80
#' nseg <- length(mp_fluss_data$walkjogrun$gtruth) # 2
#' mp <- stomp(data, window_size = w)
#' cac <- fluss_cac(mp$pi, w)
#' segments <- fluss_extract(cac, nseg, w)
#' }
fluss_extract <- function(.mpac, num_segments, exclusion_zone = NULL) {
  if (!any(class(.mpac) %in% "ArcCounts")) {
    stop("Error: First argument must be an object of class `ArcCounts`.")
  }

  if (is.null(exclusion_zone)) {
    exclusion_zone <- .mpac$ez * 10 # normally ez is 0.5, so ez here is 5
  }

  cac <- .mpac$cac # keep cac intact

  segments_positions <- vector(mode = "numeric")
  arc_counts_size <- length(cac)
  exclusion_zone <- round(.mpac$w * exclusion_zone + vars()$eps)

  for (i in 1:num_segments) {
    idx <- which.min(cac)
    if (cac[idx] >= 1) {
      break
    }
    segments_positions[i] <- idx
    cac[max(1, (idx - exclusion_zone)):min(arc_counts_size, (idx + exclusion_zone - 1))] <- Inf
  }

  .mpac$fluss <- segments_positions

  class(.mpac) <- append("Fluss", class(.mpac))

  return(.mpac)
}

#' FLUSS - Corrected Arc Counts
#'
#' Computes the arc count with edge correction (CAC).
#'
#' Original paper suggest using the classic statistical-process-control heuristic to set a threshold
#' where a semantic change may occur in CAC. This may be useful in real-time implementation as we don't
#' know in advance the number of domain changes to look for. Please check original paper (1).
#'
#' @param profile_index the profile index for arc counting.
#' @param window_size an `int`. Size of the sliding window.
#' @param exclusion_zone a `numeric`. Size of the exclusion zone, based on window size (default is `5`).
#'
#' @return Returns a companion matrix with the same size of profile index. This matrix contains the number of
#' 'arcs' crossing over each index.
#'
#' @export
#' @family fluss
#' @references * Gharghabi S, Ding Y, Yeh C-CM, Kamgar K, Ulanova L, Keogh E. Matrix Profile VIII: Domain Agnostic Online Semantic Segmentation at Superhuman Performance Levels. In: 2017 IEEE International Conference on Data Mining (ICDM). IEEE; 2017. p. 117–26.
#' @references Website: <https://sites.google.com/site/onlinesemanticsegmentation/>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#' @examples
#' # This is a fast toy example and results are useless. For a complete result, run the code inside
#' #'Not run' section below.
#' data <- mp_fluss_data$tilt_abp$data[1:1000]
#' w <- 210
#' mp <- stomp(data, window_size = w, verbose = 0)
#' cac <- fluss_cac(mp$pi, w)
#'
#' \dontrun{
#' data <- mp_fluss_data$walkjogrun$data
#' w <- mp_fluss_data$walkjogrun$window # 80
#' mp <- stomp(data, window_size = w)
#' cac <- fluss_cac(mp$pi, w)
#' }
fluss_cac <- function(.mp, exclusion_zone = NULL) {
  if (!any(class(.mp) %in% "MatrixProfile")) {
    stop("Error: First argument must be an object of class `MatrixProfile`.")
  }

  if (is.null(exclusion_zone)) {
    exclusion_zone <- .mp$ez * 10 # normally ez is 0.5, so ez here is 5
  }

  arc_counts <- vector(mode = "numeric")
  profile_index_size <- length(.mp$pi)

  nnmark <- matrix(0, profile_index_size, 1)

  for (i in 1:profile_index_size) {
    j <- .mp$pi[i]
    nnmark[min(i, j)] <- nnmark[min(i, j)] + 1
    nnmark[max(i, j)] <- nnmark[max(i, j)] - 1
  }

  arc_counts <- cumsum(nnmark)

  ideal_arc_counts <- stats::dbeta(seq(0, 1, length.out = profile_index_size), 2, 2) * profile_index_size / 3
  corrected_arc_counts <- pmin(arc_counts / ideal_arc_counts, 1)
  exclusion_zone <- round(.mp$w * exclusion_zone + vars()$eps)
  corrected_arc_counts[1:min(exclusion_zone, profile_index_size)] <- 1
  corrected_arc_counts[max((profile_index_size - exclusion_zone + 1), 1):profile_index_size] <- 1

  .mp$cac <- corrected_arc_counts

  class(.mp) <- append("ArcCounts", class(.mp))

  return(.mp)
}

#' FLUSS - Prediction score calculation
#'
#' @param gtruth an `int` or `vector` of `int` with the ground truth index of segments.
#' @param extracted an `int` or `vector` of `int` with the extracted indexes from [fluss_extract()].
#' @param data_size an `int`. Size of original input data.
#'
#' @return Returns the score of predicted semantic transitions compared with the ground truth.
#' Zero is the best, One is the worst.
#' @export
#' @family fluss
#' @references * Gharghabi S, Ding Y, Yeh C-CM, Kamgar K, Ulanova L, Keogh E. Matrix Profile VIII: Domain Agnostic Online Semantic Segmentation at Superhuman Performance Levels. In: 2017 IEEE International Conference on Data Mining (ICDM). IEEE; 2017. p. 117–26.
#' @references Website: <https://sites.google.com/site/onlinesemanticsegmentation/>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#' @examples
#' # This is a fast toy example and results are useless. For a complete result, run the code inside
#' #'Not run' section below.
#' data <- mp_fluss_data$tilt_abp$data[1:1000]
#' w <- 10
#' truth <- 400
#' mp <- stomp(data, window_size = w, verbose = 0)
#' cac <- fluss_cac(mp$pi, w)
#' segments <- fluss_extract(cac, 1, w)
#' score <- fluss_score(truth, segments, length(data))
#' \dontrun{
#' data <- mp_fluss_data$walkjogrun$data
#' w <- mp_fluss_data$walkjogrun$window # 80
#' truth <- mp_fluss_data$walkjogrun$gtruth # 3800 6800
#' nseg <- length(mp_fluss_data$walkjogrun$gtruth) # 2
#' mp <- stomp(data, window_size = w)
#' cac <- fluss_cac(mp$pi, w)
#' segments <- fluss_extract(cac, nseg, w)
#' score <- fluss_score(truth, segments, length(data))
#' }
fluss_score <- function(gtruth, extracted, data_size) {
  n <- length(gtruth)
  m <- length(extracted)
  minv <- rep(Inf, n)

  for (j in 1:n) {
    for (i in 1:m) {
      if (abs(extracted[i] - gtruth[j]) < abs(minv[j])) {
        minv[j] <- abs(extracted[i] - gtruth[j])
      }
    }
  }

  score <- sum(minv) / data_size

  return(score)
}
