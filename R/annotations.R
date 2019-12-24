#' Computes the annotation vector that favors number of zero crossing
#'
#' @param .mp a Matrix Profile object.
#' @param data a `vector` or a column `matrix` of `numeric`.
#' @param apply logical. (Default is `FALSE`). Applies the Annotation Vector over the Matrix Profile.
#'  Use with caution.
#'
#' @return Returns the input `.mp` object with an embedded annotation vector.
#' @export
#' @family Annotation vectors
#' @references * Dau HA, Keogh E. Matrix Profile V: A Generic Technique to Incorporate Domain
#'   Knowledge into Motif Discovery. In: Proceedings of the 23rd ACM SIGKDD International Conference
#'   on Knowledge Discovery and Data Mining - KDD '17. New York, New York, USA: ACM Press; 2017. p.
#'   125-34.
#' @examples
#' data <- mp_test_data$train$data[1:1000]
#' w <- 50
#' mp <- tsmp(data, window_size = w, verbose = 0)
#' av <- av_zerocrossing(mp, apply = TRUE)
av_zerocrossing <- function(.mp, data, apply = FALSE) {
  if ("Valmod" %in% class(.mp)) {
    stop("Function not implemented for objects of class `Valmod`.")
  }

  if (missing(data) && !is.null(.mp$data)) {
    data <- .mp$data[[1L]]
  }

  data <- as.matrix(data)
  data <- as.matrix(data[, 1L])
  data <- znorm(data)
  profile_size <- length(.mp$mp)
  av <- matrix(0L, profile_size, 1L)
  for (j in 1L:profile_size) {
    av[j] <- zero_crossings(data[j:(j + .mp$w - 1L), ])
  }

  av <- zero_one_norm(av)

  .mp$av <- av

  class(.mp) <- update_class(class(.mp), "AnnotationVector")

  if (apply == TRUE) {
    .mp <- av_apply(.mp)
  }

  return(.mp)
}

#' Computes the annotation vector that favors complexity
#'
#' @inheritParams av_zerocrossing
#' @param dilution_factor a `numeric`. (Default is `0`). Larger numbers means more dilution.
#'
#' @return Returns the input `.mp` object with an embedded annotation vector.
#' @export
#' @family Annotation vectors
#' @references * Dau HA, Keogh E. Matrix Profile V: A Generic Technique to Incorporate Domain
#'   Knowledge into Motif Discovery. In: Proceedings of the 23rd ACM SIGKDD International Conference
#'   on Knowledge Discovery and Data Mining - KDD '17. New York, New York, USA: ACM Press; 2017. p.
#'   125-34.
#' @examples
#' data <- mp_test_data$train$data[1:1000]
#' w <- 50
#' mp <- tsmp(data, window_size = w, verbose = 0)
#' av <- av_complexity(mp, apply = TRUE)
av_complexity <- function(.mp, data, dilution_factor = 0, apply = FALSE) {
  if ("Valmod" %in% class(.mp)) {
    stop("Function not implemented for objects of class `Valmod`.")
  }

  if (missing(data) && !is.null(.mp$data)) {
    data <- .mp$data[[1]]
  }

  data <- as.matrix(data)
  data <- as.matrix(data[, 1])
  data <- znorm(data)
  profile_size <- length(.mp$mp)
  av <- matrix(0, profile_size, 1)

  for (j in 1:profile_size) {
    av[j] <- complexity(data[j:(j + .mp$w - 1)])
  }

  av <- zero_one_norm(av) # zero-one normalize the av

  # Select dilution factor, 0 is no dilution,
  # larger numbers are more dilution
  av <- av + dilution_factor
  av <- av / (dilution_factor + 1L)

  .mp$av <- av

  class(.mp) <- update_class(class(.mp), "AnnotationVector")

  if (apply == TRUE) {
    .mp <- av_apply(.mp)
  }

  return(.mp)
}


#' Computes the annotation vector that suppresses motion artifacts
#'
#' @inheritParams av_zerocrossing
#'
#' @return Returns the input `.mp` object with an embedded annotation vector.
#' @export
#' @family Annotation vectors
#' @references * Dau HA, Keogh E. Matrix Profile V: A Generic Technique to Incorporate Domain
#'   Knowledge into Motif Discovery. In: Proceedings of the 23rd ACM SIGKDD International Conference
#'   on Knowledge Discovery and Data Mining - KDD '17. New York, New York, USA: ACM Press; 2017. p.
#'   125-34.
#' @examples
#' data <- mp_test_data$train$data[1:1000]
#' w <- 50
#' mp <- tsmp(data, window_size = w, verbose = 0)
#' av <- av_motion_artifact(mp, apply = TRUE)
av_motion_artifact <- function(.mp, data, apply = FALSE) {
  if ("Valmod" %in% class(.mp)) {
    stop("Function not implemented for objects of class `Valmod`.")
  }

  if (missing(data) && !is.null(.mp$data)) {
    data <- .mp$data[[1]]
  }

  data <- as.matrix(data)
  data <- as.matrix(data[, 1])
  data <- znorm(data)
  profile_size <- length(.mp$mp)
  av <- matrix(0L, profile_size, 1)

  for (i in 1:profile_size) {
    s <- data[i:(i + .mp$w - 1L), ]
    av[i] <- stats::sd(s)
  }

  cav <- av
  mu <- mean(av)

  cav[av >= mu] <- 0
  cav[av < mu] <- 1

  .mp$av <- cav

  class(.mp) <- update_class(class(.mp), "AnnotationVector")

  if (apply == TRUE) {
    .mp <- av_apply(.mp)
  }

  return(.mp)
}

#' Computes the annotation vector that suppresses stop-word motifs
#'
#' @details
#' The function is intended to be generic. However, its parameters (`stop_word_loc`,
#' `exclusion_zone` and `threshold`) are highly dataset dependent.
#'
#' @inheritParams av_zerocrossing
#' @param stop_word_loc an `int`. The index of stop word location.
#' @param exclusion_zone a `numeric`. Size of the exclusion zone, based on window_size (default is
#'   `NULL`). See details.
#' @param threshold a `numeric`. (default is `0.1`).
#'
#' @return Returns the input `.mp` object with an embedded annotation vector.
#' @export
#' @family Annotation vectors
#' @references * Dau HA, Keogh E. Matrix Profile V: A Generic Technique to Incorporate Domain
#'   Knowledge into Motif Discovery. In: Proceedings of the 23rd ACM SIGKDD International Conference
#'   on Knowledge Discovery and Data Mining - KDD '17. New York, New York, USA: ACM Press; 2017. p.
#'   125-34.
#' @examples
#' data <- mp_test_data$train$data[1:1000]
#' w <- 50
#' mp <- tsmp(data, window_size = w, verbose = 0)
#' av <- av_stop_word(mp, stop_word_loc = 150, apply = TRUE)
av_stop_word <- function(.mp, data, stop_word_loc, exclusion_zone = NULL, threshold = 0.1, apply = FALSE) {
  if ("Valmod" %in% class(.mp)) {
    stop("Function not implemented for objects of class `Valmod`.")
  }

  if (missing(data) && !is.null(.mp$data)) {
    data <- .mp$data[[1L]]
  }

  if (is.null(exclusion_zone)) {
    exclusion_zone <- .mp$ez
  }

  data <- as.matrix(data)
  data <- as.matrix(data[, 1])
  data <- znorm(data)
  stop_word <- data[stop_word_loc:(stop_word_loc + .mp$w - 1L), ]

  profile_size <- length(.mp$mp)

  av <- matrix(0L, profile_size, 1L)

  for (i in 1L:profile_size) {
    s <- data[i:(i + .mp$w - 1L), ]
    av[i, ] <- diff2(s, stop_word)
  }

  av <- zero_one_norm(av)

  index <- which(av <= threshold)

  for (i in seq_len(length(index))) {
    if (index[i] < exclusion_zone) {
      av[(index[i] - index[i] + 1):min((index[i] + exclusion_zone - 1L), profile_size), ] <- 0L
    } else {
      av[(index[i] - exclusion_zone + 1L):min((index[i] + exclusion_zone - 1L), profile_size), ] <- 0L
    }
  }

  .mp$av <- av

  class(.mp) <- update_class(class(.mp), "AnnotationVector")

  if (apply == TRUE) {
    .mp <- av_apply(.mp)
  }

  return(.mp)
}

#' Computes the annotation vector that suppresses hard-limited artifacts
#'
#' @inheritParams av_zerocrossing
#'
#' @return Returns the input `.mp` object with an embedded annotation vector.
#' @export
#' @family Annotation vectors
#' @references * Dau HA, Keogh E. Matrix Profile V: A Generic Technique to Incorporate Domain
#'   Knowledge into Motif Discovery. In: Proceedings of the 23rd ACM SIGKDD International Conference
#'   on Knowledge Discovery and Data Mining - KDD '17. New York, New York, USA: ACM Press; 2017. p.
#'   125-34.
#' @examples
#' data <- mp_test_data$train$data[1:1000]
#' w <- 50
#' mp <- tsmp(data, window_size = w, verbose = 0)
#' av <- av_hardlimit_artifact(mp, apply = TRUE)
av_hardlimit_artifact <- function(.mp, data, apply = FALSE) {
  if ("Valmod" %in% class(.mp)) {
    stop("Function not implemented for objects of class `Valmod`.")
  }

  if (missing(data) && !is.null(.mp$data)) {
    data <- .mp$data[[1L]]
  }

  data <- as.matrix(data)
  data <- as.matrix(data[, 1L])
  data <- znorm(data)
  max <- max(data)
  min <- min(data)

  profile_size <- length(.mp$mp)
  av <- matrix(0, profile_size, 1L)

  for (i in 1L:profile_size) {
    s <- data[i:(i + .mp$w - 1L), ]
    av[i, ] <- length(s[s == max | s == min])
  }

  av <- zero_one_norm(av) # zero-one normalize the av
  av <- 1L - av

  .mp$av <- av

  class(.mp) <- update_class(class(.mp), "AnnotationVector")

  if (apply == TRUE) {
    .mp <- av_apply(.mp)
  }

  return(.mp)
}

#' Corrects the matrix profile using an annotation vector
#'
#' This function overwrites the current Matrix Profile using the Annotation Vector. Use with caution.
#'
#' @param .mp A Matrix Profile with an Annotation Vector.
#'
#' @return Returns the input `.mp` object corrected by the embedded annotation vector.
#' @export
#' @family Annotation vectors
#' @references * Dau HA, Keogh E. Matrix Profile V: A Generic Technique to Incorporate Domain
#'   Knowledge into Motif Discovery. In: Proceedings of the 23rd ACM SIGKDD International Conference
#'   on Knowledge Discovery and Data Mining - KDD '17. New York, New York, USA: ACM Press; 2017. p.
#'   125-34.
#' @examples
#' data <- mp_test_data$train$data[1:1000]
#' w <- 50
#' mp <- tsmp(data, window_size = w, verbose = 0)
#' mp <- av_complexity(mp)
#' av <- av_apply(mp)
av_apply <- function(.mp) {
  if (!("MatrixProfile" %in% class(.mp))) {
    stop("First argument must be an object of class `MatrixProfile`.")
  }

  if (!("AnnotationVector" %in% class(.mp))) {
    stop("First argument must be an object of class `AnnotationVector`.")
  }

  if (!is.null(attr(.mp, "annotated"))) {
    stop("This Matrix Profile has already been annotated.")
  }

  if (sys.parent() == 0) {
    warning("Warning: This function overwrites the current Matrix Profile.")
  }

  .mp$mp <- .mp$mp + (1 - .mp$av) * max(.mp$mp)

  attr(.mp, "annotated") <- TRUE

  return(.mp)
}
