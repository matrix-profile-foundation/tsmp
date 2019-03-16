#' Fast Low-cost Unipotent Semantic Segmentation (FLUSS)
#'
#' FLUSS is a Domain Agnostic Online Semantic Segmentation that uses the assumption that when few
#' arc are crossing a given index point, means that there is a high probability of semantic change.
#' This function is a wrap to [fluss_cac()] and [fluss_extract()].
#'
#' @param .mp a TSMP object of class `MatrixProfile`.
#' @param num_segments an `int`. Number of segments to extract. Based on domain knowledge.
#' @param exclusion_zone if a `number` will be used instead of embedded value. (Default is `NULL`).
#'
#' @return Returns the input `.mp` object new names: `cac`, corrected arc count and `fluss` with
#' the location of semantic changes.
#' @export
#' @family Semantic Segmentations
#' @references * Gharghabi S, Ding Y, Yeh C-CM, Kamgar K, Ulanova L, Keogh E. Matrix Profile VIII:
#'   Domain Agnostic Online Semantic Segmentation at Superhuman Performance Levels. In: 2017 IEEE
#'   International Conference on Data Mining (ICDM). IEEE; 2017. p. 117–26.
#' @references Website: <https://sites.google.com/site/onlinesemanticsegmentation/>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#' @examples
#' data <- mp_fluss_data$tilt_abp$data[1:1000]
#' w <- 10
#' mp <- tsmp(data, window_size = w, verbose = 0)
#' mp <- fluss(mp, 2)
fluss <- function(.mp, num_segments, exclusion_zone = NULL) {
  if ("Valmod" %in% class(.mp)) {
    stop("Function not implemented for objects of class `Valmod`.")
  }

  fluss_extract(fluss_cac(.mp, exclusion_zone), num_segments = num_segments, exclusion_zone = exclusion_zone)
}

#' FLUSS - Extract Segments
#'
#' Extract candidate points of semantic changes.
#'
#' @param .mpac a TSMP object of class `ArcCount`.
#' @param num_segments an `int`. Number of segments to extract. Based on domain knowledge.
#' @param exclusion_zone if a `number` will be used instead of embedded value. (Default is `NULL`).
#'
#' @return Returns the input `.mp` object a new name `fluss` with the location of semantic changes.
#' @export
#' @family Semantic Segmentations
#' @references * Gharghabi S, Ding Y, Yeh C-CM, Kamgar K, Ulanova L, Keogh E. Matrix Profile VIII:
#'   Domain Agnostic Online Semantic Segmentation at Superhuman Performance Levels. In: 2017 IEEE
#'   International Conference on Data Mining (ICDM). IEEE; 2017. p. 117–26.
#' @references Website: <https://sites.google.com/site/onlinesemanticsegmentation/>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#' @examples
#' data <- mp_fluss_data$tilt_abp$data[1:1000]
#' w <- 10
#' mp <- tsmp(data, window_size = w, verbose = 0)
#' mp <- fluss_cac(mp)
#' mp <- fluss_extract(mp, 2)
fluss_extract <- function(.mpac, num_segments, exclusion_zone = NULL) {
  if (!any(class(.mpac) %in% "ArcCount")) {
    stop("First argument must be an object of class `ArcCount`.")
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

  class(.mpac) <- update_class(class(.mpac), "Fluss")

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
#' @param .mp a TSMP object of class `MatrixProfile`.
#' @param exclusion_zone if a `number` will be used instead of embedded value. (Default is `NULL`).
#'
#' @return Returns the input `.mp` object a new name `cac` with the corrected arc count.
#'
#' @export
#' @family Semantic Segmentations
#' @references * Gharghabi S, Ding Y, Yeh C-CM, Kamgar K, Ulanova L, Keogh E. Matrix Profile VIII: Domain Agnostic Online Semantic Segmentation at Superhuman Performance Levels. In: 2017 IEEE International Conference on Data Mining (ICDM). IEEE; 2017. p. 117–26.
#' @references Website: <https://sites.google.com/site/onlinesemanticsegmentation/>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#' @examples
#' data <- mp_fluss_data$tilt_abp$data[1:1000]
#' w <- 10
#' mp <- tsmp(data, window_size = w, verbose = 0)
#' mp <- fluss_cac(mp)
fluss_cac <- function(.mp, exclusion_zone = NULL) {
  if (!("MatrixProfile" %in% class(.mp))) {
    stop("First argument must be an object of class `MatrixProfile`.")
  }

  if ("Valmod" %in% class(.mp)) {
    stop("Function not implemented for objects of class `Valmod`.")
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

  class(.mp) <- update_class(class(.mp), "ArcCount")

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
#' @family Semantic Segmentations
#' @references * Gharghabi S, Ding Y, Yeh C-CM, Kamgar K, Ulanova L, Keogh E. Matrix Profile VIII: Domain Agnostic Online Semantic Segmentation at Superhuman Performance Levels. In: 2017 IEEE International Conference on Data Mining (ICDM). IEEE; 2017. p. 117–26.
#' @references Website: <https://sites.google.com/site/onlinesemanticsegmentation/>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#' @examples
#' data <- mp_fluss_data$tilt_abp$data[1:1000]
#' w <- 10
#' truth <- c(945, 875)
#' mp <- tsmp(data, window_size = w, verbose = 0)
#' mp <- fluss_cac(mp)
#' mp <- fluss_extract(mp, 2)
#' score <- fluss_score(truth, mp$fluss, length(data))
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
