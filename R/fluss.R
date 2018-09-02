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
#' @param window.size an `int`. Size of the sliding window.
#' @param num.segments an `int`. Number of segments to extract. Based on domain knowledge.
#' @param exclusion.zone a `numeric`. Size of the exclusion zone, based on window size (default is
#'   `5`).
#' @param gtruth an `int` or `vector` of `int` with the ground truth index of segments. (Default is
#'   `NULL`).
#' @param profile.index a pre-computed profile index. (Default is `NULL`).
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
#' data <- fluss_data$tilt.abp$data[1:1000]
#' truth <- 400
#' w <- 10
#' segments <- fluss(data, w, 1, gtruth = truth, verbose = 0)
#' \dontrun{
#' data <- fluss_data$walkjogrun$data
#' w <- fluss_data$walkjogrun$window # 80
#' truth <- fluss_data$walkjogrun$gtruth # 3800 6800
#' nseg <- length(fluss_data$walkjogrun$gtruth) # 2
#' segments <- fluss(data, w, nseg, gtruth = truth)
#' }
fluss <- function(data, window.size, num.segments, exclusion.zone = 5, gtruth = NULL, profile.index = NULL, verbose = 2) {

  ## Input validatin
  if (is.matrix(data) || is.data.frame(data)) {
    if (is.data.frame(data)) {
      data <- as.matrix(data)
    } # just to be uniform
    if (ncol(data) > nrow(data)) {
      data <- t(data)
    }
    data.size <- nrow(data)
  } else if (is.vector(data)) {
    data.size <- length(data)
    data <- as.matrix(data) # just to be uniform
  } else {
    stop("Error: Unknown type of data. Must be: matrix, data.frame or vector.", call. = FALSE)
  }

  profile <- NULL
  if (is.null(profile.index)) {
    profile <- stomp.par(data, window.size = window.size, verbose = verbose)
    profile.index <- profile$pi
  }

  cac <- fluss.cac(profile.index, window.size, exclusion.zone)

  segments <- fluss.extract(cac, num.segments, window.size, exclusion.zone)

  if (!is.null(gtruth)) {
    score <- fluss.score(gtruth, segments, data.size)
    return(list(segments = segments, mp = profile$mp, pi = profile.index, cac = cac, score = score))
  }
  else {
    return(list(segments = segments, mp = profile$mp, pi = profile.index, cac = cac))
  }
}

#' FLUSS - Extract Segments
#'
#' Extract candidate points of semantic changes.
#'
#' @param arc.counts a `matrix` with the corrected arc counts from [fluss.cac()].
#' @param num.segments an `int`. Number of segments to extract. Based on domain knowledge.
#' @param window.size an `int`. Size of the sliding window.
#' @param exclusion.zone a `numeric`. Size of the exclusion zone, based on window size (default is
#'   `5`).
#'
#' @return Returns an `int` or a `vector` of `int` with the location of predicted semantic changes.
#'   The number of locations is not greater than `num.segments`.
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
#' data <- fluss_data$tilt.abp$data[1:1000]
#' w <- 210
#' mp <- stomp(data, window.size = w, verbose = 0)
#' cac <- fluss.cac(mp$pi, w)
#' segments <- fluss.extract(cac, 1, w)
#' \dontrun{
#' data <- fluss_data$walkjogrun$data
#' w <- fluss_data$walkjogrun$window # 80
#' nseg <- length(fluss_data$walkjogrun$gtruth) # 2
#' mp <- stomp(data, window.size = w)
#' cac <- fluss.cac(mp$pi, w)
#' segments <- fluss.extract(cac, nseg, w)
#' }
fluss.extract <- function(arc.counts, num.segments, window.size, exclusion.zone = 5) {
  segments.positions <- vector(mode = "numeric")
  arc.counts.size <- length(arc.counts)
  exclusion.zone <- floor(window.size * exclusion.zone)

  for (i in 1:num.segments) {
    idx <- which.min(arc.counts)
    if (arc.counts[idx] >= 1) {
      break
    }
    segments.positions[i] <- idx
    arc.counts[max(1, (idx - exclusion.zone)):min(arc.counts.size, (idx + exclusion.zone - 1))] <- Inf
  }

  return(segments.positions)
}

#' FLUSS - Corrected Arc Counts
#'
#' Computes the arc count with edge correction (CAC).
#'
#' Original paper suggest using the classic statistical-process-control heuristic to set a threshold
#' where a semantic change may occur in CAC. This may be useful in real-time implementation as we don't
#' know in advance the number of domain changes to look for. Please check original paper (1).
#'
#' @param profile.index the profile index for arc counting.
#' @param window.size an `int`. Size of the sliding window.
#' @param exclusion.zone a `numeric`. Size of the exclusion zone, based on window size (default is `5`).
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
#' data <- fluss_data$tilt.abp$data[1:1000]
#' w <- 210
#' mp <- stomp(data, window.size = w, verbose = 0)
#' cac <- fluss.cac(mp$pi, w)
#'
#' \dontrun{
#' data <- fluss_data$walkjogrun$data
#' w <- fluss_data$walkjogrun$window # 80
#' mp <- stomp(data, window.size = w)
#' cac <- fluss.cac(mp$pi, w)
#' }
fluss.cac <- function(profile.index, window.size, exclusion.zone = 5) {
  arc.counts <- vector(mode = "numeric")
  profile.index.size <- length(profile.index)

  nnmark <- matrix(0, profile.index.size, 1)

  for (i in 1:profile.index.size) {
    j <- profile.index[i]
    nnmark[min(i, j)] <- nnmark[min(i, j)] + 1
    nnmark[max(i, j)] <- nnmark[max(i, j)] - 1
  }

  arc.counts <- cumsum(nnmark)

  ideal.arc.counts <- stats::dbeta(seq(0, 1, length.out = profile.index.size), 2, 2) * profile.index.size / 3
  corrected.arc.counts <- pmin(arc.counts / ideal.arc.counts, 1)
  exclusion.zone <- floor(window.size * exclusion.zone)
  corrected.arc.counts[1:min(exclusion.zone, profile.index.size)] <- 1
  corrected.arc.counts[max((profile.index.size - exclusion.zone + 1), 1):profile.index.size] <- 1

  return(corrected.arc.counts)
}

#' FLUSS - Prediction score calculation
#'
#' @param gtruth an `int` or `vector` of `int` with the ground truth index of segments.
#' @param extracted an `int` or `vector` of `int` with the extracted indexes from [fluss.extract()].
#' @param data.size an `int`. Size of original input data.
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
#' data <- fluss_data$tilt.abp$data[1:1000]
#' w <- 10
#' truth <- 400
#' mp <- stomp(data, window.size = w, verbose = 0)
#' cac <- fluss.cac(mp$pi, w)
#' segments <- fluss.extract(cac, 1, w)
#' score <- fluss.score(truth, segments, length(data))
#' \dontrun{
#' data <- fluss_data$walkjogrun$data
#' w <- fluss_data$walkjogrun$window # 80
#' truth <- fluss_data$walkjogrun$gtruth # 3800 6800
#' nseg <- length(fluss_data$walkjogrun$gtruth) # 2
#' mp <- stomp(data, window.size = w)
#' cac <- fluss.cac(mp$pi, w)
#' segments <- fluss.extract(cac, nseg, w)
#' score <- fluss.score(truth, segments, length(data))
#' }
fluss.score <- function(gtruth, extracted, data.size) {
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

  score <- sum(minv) / data.size

  return(score)
}
