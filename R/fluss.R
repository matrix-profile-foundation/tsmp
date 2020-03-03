#' Fast Low-cost Unipotent Semantic Segmentation (FLUSS)
#'
#' FLUSS is a Domain Agnostic Online Semantic Segmentation that uses the assumption that when few
#' arc are crossing a given index point, means that there is a high probability of semantic change.
#' This function is a wrap to [fluss_cac()] and [fluss_extract()].
#'
#' @param .mp a `MatrixProfile` object.
#' @param num_segments an `int`. Number of segments to extract. Based on domain knowledge.
#' @param exclusion_zone if a `number` will be used instead of embedded value. (Default is `NULL`).
#'
#' @return Returns the input `.mp` object new names: `cac`, corrected arc count and `fluss` with
#' the location of semantic changes.
#' @export
#' @family Semantic Segmentations
#' @references * Gharghabi S, Ding Y, Yeh C-CM, Kamgar K, Ulanova L, Keogh E. Matrix Profile VIII:
#'   Domain Agnostic Online Semantic Segmentation at Superhuman Performance Levels. In: 2017 IEEE
#'   International Conference on Data Mining (ICDM). IEEE; 2017. p. 117-26.
#' @references Website: <https://sites.google.com/site/onlinesemanticsegmentation/>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#' @examples
#' data <- mp_fluss_data$tilt_abp$data[1:1000]
#' w <- 10
#' mp <- tsmp(data, window_size = w, verbose = 0)
#' mp <- fluss(mp, 2)
fluss <- function(.mp, num_segments = 1, exclusion_zone = NULL) {
  if ("Valmod" %in% class(.mp)) {
    stop("Function not implemented for objects of class `Valmod`.")
  }

  fluss_extract(fluss_cac(.mp, exclusion_zone), num_segments = num_segments, exclusion_zone = exclusion_zone)
}

#' Fast Low-cost Online Semantic Segmentation (FLOSS)
#'
#' @param .mp a `MatrixProfile` object.
#' @param new_data a `matrix`or `vector` of new observations.
#' @param data_window an `int`. Sets the size of the buffer used to keep track of semantic changes.
#' @param threshold a `number`. (Default is `1`). Set the maximum value for evaluating semantic changes.
#' This is data specific. It is advised to check what is 'normal' for your data.
#' @param exclusion_zone if a `number` will be used instead of embedded value. (Default is `NULL`).
#' @param chunk_size an `int` . (Default is `NULL`). Set the size of new data that will be added to
#' Floss in each iteration if `new_data` is large. If `NULL`, the size will be 50. This is not needed
#' if `new_data` is small, like 1 observation.
#' @param keep_cac a `logical`. (Default is `TRUE`). If set to `FALSE`, the `cac_final` will contain
#' only values within `data_window`
#'
#' @return Returns the input `.mp` object new names: `cac` the corrected arc count, `cac_final`the
#' combination of `cac` after repeated calls of `floss()`, `floss` with the location of semantic
#' changes and `floss_vals` with the normalized arc count value of the semantic change positions.
#' @export
#' @family Semantic Segmentations
#' @references * Gharghabi S, Ding Y, Yeh C-CM, Kamgar K, Ulanova L, Keogh E. Matrix Profile VIII:
#'   Domain Agnostic Online Semantic Segmentation at Superhuman Performance Levels. In: 2017 IEEE
#'   International Conference on Data Mining (ICDM). IEEE; 2017. p. 117-26.
#' @references Website: <https://sites.google.com/site/onlinesemanticsegmentation/>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#' @examples
#' data <- mp_fluss_data$tilt_abp$data[1:1000]
#' new_data <- mp_fluss_data$tilt_abp$data[1001:1010]
#' new_data2 <- mp_fluss_data$tilt_abp$data[1011:1020]
#' w <- 80
#' mp <- tsmp(data, window_size = w, verbose = 0)
#' data_window <- 1000
#' mp <- floss(mp, new_data, data_window)
#' mp <- floss(mp, new_data2, data_window)
floss <- function(.mp, new_data, data_window, threshold = 1, exclusion_zone = NULL, chunk_size = NULL, keep_cac = TRUE) {
  if (missing(data_window)) {
    stop("argument 'data_window' is missing, looking for fluss() instead?")
  }

  if ("Valmod" %in% class(.mp)) {
    stop("Function not implemented for objects of class `Valmod`.")
  }

  mp_offset <- ifelse(!is.null(attr(.mp, "offset")), attr(.mp, "offset"), 0)
  data_size <- nrow(.mp$data[[1]])
  new_data_size <- length(new_data)

  if (data_size < data_window) {
    if ((data_size + new_data_size) <= data_window) {
      .mp <- stompi_update(.mp, new_data)
      .mp$cac_final <- NULL
      return(.mp)
    } else {
      .mp <- stompi_update(.mp, utils::head(new_data, data_window - data_size))
      new_data <- new_data[(data_window - data_size + 1):new_data_size]
    }
  }

  if (is.null(chunk_size)) {
    chunk_size <- min(floor(nrow(.mp$data[[1]]) / 2), 50)
    chunk_size <- min(chunk_size, floor(data_window / 2))
  }

  if (data_size > data_window) {
    if (mp_offset > 0) {
      attr(.mp, "new_data") <- 0
      .mp <- floss_cac(.mp, data_size, 0)
    }
    else {
      .mp <- fluss_cac(.mp, 0)
    }

    .mp$cac_final <- c(
      rep(NA, chunk_size + mp_offset - chunk_size),
      utils::head(.mp$cac, -(round(data_window * (1 - vars()$kmode) - (1 - vars()$kmode) * .mp$w) - 0.5 * chunk_size + ifelse(mp_offset > 0, 1, 2)))
    )

    na_head <- round(vars()$kmode * data_window + (0.5 * chunk_size - vars()$kmode * .mp$w)) + mp_offset

    .mp$cac_final[1:na_head] <- NA
  }

  num_chunks <- floor(length(new_data) / chunk_size)
  last_chunk <- length(new_data) - num_chunks * chunk_size

  end_idx <- 0
  for (i in seq_len(num_chunks)) {
    st_idx <- chunk_size * i - chunk_size + 1
    end_idx <- st_idx + chunk_size - 1
    .mp <- floss_cac(stompi_update(.mp, new_data[st_idx:end_idx], data_window), data_window, exclusion_zone)
  }

  if (last_chunk > 0) {
    s_idx <- end_idx + 1
    e_idx <- s_idx + last_chunk - 1
    .mp <- floss_cac(stompi_update(.mp, new_data[s_idx:e_idx], data_window), data_window, exclusion_zone)
  }

  res <- floss_extract(.mp, threshold)

  if (!keep_cac) {
    res$cac_final <- utils::tail(res$cac_final, -length(new_data))
  }

  return(res)
}

#' FLOSS - Extract Segments
#'
#' Extract candidate points of semantic changes.
#'
#' @param .mpac a TSMP object of class `ArcCount`.
#' @param exclusion_zone if a `number` will be used instead of embedded value. (Default is `NULL`).
#' @param threshold a `number`. (Default is `1`). Set the maximum value for evaluating semantic changes.
#' This is data specific. It is advised to check what is 'normal' for your data.
#'
#' @return Returns the input `.mp` object a new name `floss` with the location of semantic
#' changes and `floss_vals` with the normalized arc count value of the semantic change positions.
#' @export
#' @family Semantic Segmentations
#' @references * Gharghabi S, Ding Y, Yeh C-CM, Kamgar K, Ulanova L, Keogh E. Matrix Profile VIII:
#'   Domain Agnostic Online Semantic Segmentation at Superhuman Performance Levels. In: 2017 IEEE
#'   International Conference on Data Mining (ICDM). IEEE; 2017. p. 117-26.
#' @references Website: <https://sites.google.com/site/onlinesemanticsegmentation/>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#' @examples
#' data <- mp_fluss_data$tilt_abp$data[1:1000]
#' w <- 10
#' mp <- tsmp(data, window_size = w, verbose = 0)
#' mp <- fluss_cac(mp)
#' mp <- fluss_extract(mp, 2)
floss_extract <- function(.mpac, threshold = 1, exclusion_zone = NULL) {
  if (!any(class(.mpac) %in% "ArcCount")) {
    stop("First argument must be an object of class `ArcCount`.")
  }

  if (is.null(.mpac$cac_final)) {
    stop("There is no real-time information to extract. Looking for fluss_extract() instead?")
  }

  if (is.null(exclusion_zone)) {
    exclusion_zone <- .mpac$ez * 10 # normally ez is 0.5, so ez here is 5
  }

  offset <- ifelse(is.null(attr(.mpac, "offset")), 0, attr(.mpac, "offset"))
  cac_fin_len <- length(.mpac$cac_final)
  mp_len <- nrow(.mpac$mp)
  new_data <- ifelse(is.null(attr(.mpac, "new_data")), 0, attr(.mpac, "new_data"))

  if (offset == 0 || cac_fin_len == floor((mp_len * vars()$kmode + new_data * 1.5))) {
    cac <- utils::tail(.mpac$cac_final, -new_data)
  } else {
    cac <- utils::tail(.mpac$cac_final, -offset)
  }

  cac[cac > threshold] <- NA

  segments <- .mpac$floss
  seg_vals <- .mpac$floss_vals
  exclusion_zone <- round(.mpac$w * exclusion_zone + vars()$eps)

  idx <- which.min(cac)

  if (length(idx) > 0) {
    val <- cac[idx]

    real_idx <- idx + offset
    seg_len <- length(segments)

    if (seg_len > 0) {
      last_idx <- segments[seg_len]
      last_val <- seg_vals[seg_len]

      if (real_idx > last_idx) {
        if (real_idx < (last_idx + exclusion_zone)) {
          # update split
          if (val < last_val) {
            segments[seg_len] <- real_idx
            seg_vals[seg_len] <- val
          }
        } else {
          # new split
          segments <- c(segments, real_idx)
          seg_vals <- c(seg_vals, val)
        }
      }
    } else {
      segments <- real_idx
      seg_vals <- val
    }
  }

  .mpac$floss <- segments
  .mpac$floss_vals <- seg_vals

  class(.mpac) <- update_class(class(.mpac), "Floss")

  return(.mpac)
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
#'   International Conference on Data Mining (ICDM). IEEE; 2017. p. 117-26.
#' @references Website: <https://sites.google.com/site/onlinesemanticsegmentation/>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#' @examples
#' data <- mp_fluss_data$tilt_abp$data[1:1000]
#' w <- 10
#' mp <- tsmp(data, window_size = w, verbose = 0)
#' mp <- fluss_cac(mp)
#' mp <- fluss_extract(mp, 2)
fluss_extract <- function(.mpac, num_segments = 1, exclusion_zone = NULL) {
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
#' @param .mp a `MatrixProfile` object.
#' @param exclusion_zone if a `number` will be used instead of embedded value. (Default is `NULL`).
#'
#' @return Returns the input `.mp` object a new name `cac` with the corrected arc count.
#'
#' @export
#' @family Semantic Segmentations
#' @references * Gharghabi S, Ding Y, Yeh C-CM, Kamgar K, Ulanova L, Keogh E. Matrix Profile VIII: Domain Agnostic Online Semantic Segmentation at Superhuman Performance Levels. In: 2017 IEEE International Conference on Data Mining (ICDM). IEEE; 2017. p. 117-26.
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

    if (j < 0 || j > profile_index_size) {
      next
    }

    nnmark[min(i, j)] <- nnmark[min(i, j)] + 1
    nnmark[max(i, j)] <- nnmark[max(i, j)] - 1
  }

  arc_counts <- cumsum(nnmark)

  x <- seq(0, 1, length.out = profile_index_size)

  if (is.null(attr(.mp, "subsetting"))) {
    ideal_arc_counts <- stats::dbeta(x, 2, 2) * profile_index_size / 3
  } else {
    ideal_arc_counts <- stats::dbeta(x, 2.1, 2.1) * profile_index_size / 3
  }

  corrected_arc_counts <- pmin(arc_counts / ideal_arc_counts, 1)
  exclusion_zone <- round(.mp$w * exclusion_zone + vars()$eps)
  corrected_arc_counts[1:min(exclusion_zone, profile_index_size)] <- 1
  corrected_arc_counts[max((profile_index_size - exclusion_zone + 1), 1):profile_index_size] <- 1

  .mp$cac <- corrected_arc_counts

  class(.mp) <- update_class(class(.mp), "ArcCount")

  return(.mp)
}

#' FLOSS - Corrected Arc Counts
#'
#' Computes the arc count with edge and 'online' correction (CAC).
#'
#' Original paper suggest using the classic statistical-process-control heuristic to set a threshold
#' where a semantic change may occur in CAC. This may be useful in real-time implementation as we don't
#' know in advance the number of domain changes to look for. Please check original paper (1).
#'
#' @param .mp a `MatrixProfile` object.
#' @param data_window an `int`. Sets the size of the buffer used to keep track of semantic changes.
#' @param exclusion_zone if a `number` will be used instead of embedded value. (Default is `NULL`).
#'
#' @return Returns the input `.mp` object a new name `cac` with the corrected arc count and `cac_final`
#' the combination of `cac` after repeated calls of `floss()`.
#' @export
#' @family Semantic Segmentations
#' @references * Gharghabi S, Ding Y, Yeh C-CM, Kamgar K, Ulanova L, Keogh E. Matrix Profile VIII: Domain Agnostic Online Semantic Segmentation at Superhuman Performance Levels. In: 2017 IEEE International Conference on Data Mining (ICDM). IEEE; 2017. p. 117-26.
#' @references Website: <https://sites.google.com/site/onlinesemanticsegmentation/>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#' @examples
#' data <- mp_fluss_data$tilt_abp$data[1:1000]
#' new_data <- mp_fluss_data$tilt_abp$data[1001:1010]
#' w <- 10
#' mp <- tsmp(data, window_size = w, verbose = 0)
#' data_window <- 1000
#' mp <- stompi_update(mp, new_data, data_window)
#' mp <- floss_cac(mp, data_window)
floss_cac <- function(.mp, data_window, exclusion_zone = NULL) {
  if (!("MatrixProfile" %in% class(.mp))) {
    stop("First argument must be an object of class `MatrixProfile`.")
  }

  if ("Valmod" %in% class(.mp)) {
    stop("Function not implemented for objects of class `Valmod`.")
  }

  if (data_window <= .mp$w) {
    stop("data_window must be larger than matrix profile's window_size: ", .mp$w)
  }

  profile_size <- nrow(.mp$mp)
  cac_size <- data_window - .mp$w + 1
  start_idx <- profile_size - cac_size + 1
  new_data_size <- attr(.mp, "new_data")
  new_data_size <- ifelse(is.null(new_data_size), 0, new_data_size)
  mp_offset <- attr(.mp, "offset")
  mp_offset <- ifelse(is.null(mp_offset), 0, mp_offset)

  exclusion_zone <- round(.mp$w * .mp$ez + vars()$eps)
  end_idx <- profile_size - exclusion_zone - 1
  pi <- .mp$pi[start_idx:end_idx]

  nnmark <- matrix(0, cac_size, 1)

  for (i in seq_along(pi)) {
    j <- pi[i]
    if (j < 0 || j > cac_size) {
      next
    }

    nnmark[min(i, j)] <- nnmark[min(i, j)] + 1
    nnmark[max(i, j)] <- nnmark[max(i, j)] - 1
  }

  arc_counts <- cumsum(nnmark)
  x <- seq(0, 1, length.out = cac_size)

  if (mp_offset > 0) {
    mode <- vars()$kmode
    a <- 1.939274 # If you change this, change vars()$kmode
    b <- 1.698150 # If you change this, change vars()$kmode
    ideal_arc_counts <- a * b * x^(a - 1) * (1 - x^a)^(b - 1) * cac_size / 4.035477 # kumaraswamy distribution
  } else {
    mode <- 0.5
    ideal_arc_counts <- stats::dbeta(x, 2, 2) * cac_size / 3
  }

  corrected_arc_counts <- pmin(arc_counts / ideal_arc_counts, 1)
  corrected_arc_counts[1:min(exclusion_zone, cac_size)] <- 1
  corrected_arc_counts[corrected_arc_counts < 0 | is.na(corrected_arc_counts)] <- 1
  mid_idx <- round(cac_size * mode) - floor(new_data_size / 2) # same as which.max(ideal_arc_counts)

  if (is.null(.mp$cac_final)) {
    # TODO: get a more reliable value than attr(.mp, "origin")$data_size
    .mp$cac_final <- rep(NA, round(data_window * (mode - 1) - new_data_size / 2 +
      max(nrow(.mp$data[[1]]), attr(.mp, "origin")$data_size) -
      (.mp$w * mode)) + mp_offset)
  }

  .mp$cac_final <- c(.mp$cac_final, corrected_arc_counts[mid_idx:(mid_idx + new_data_size - 1)])
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
#' @references * Gharghabi S, Ding Y, Yeh C-CM, Kamgar K, Ulanova L, Keogh E. Matrix Profile VIII: Domain Agnostic Online Semantic Segmentation at Superhuman Performance Levels. In: 2017 IEEE International Conference on Data Mining (ICDM). IEEE; 2017. p. 117-26.
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
