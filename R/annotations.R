#' Computes the annotation vector that favors complexity
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#' @param window_size an `int`. Size of the sliding window.
#' @param dilution_factor a `numeric`. (Default is `0`). Larger numbers means more dilution.
#'
#' @return Returns the annotation vector for matrix profile correction.
#' @export
#' @references * Dau HA, Keogh E. Matrix Profile V: A Generic Technique to Incorporate Domain
#'   Knowledge into Motif Discovery. In: Proceedings of the 23rd ACM SIGKDD International Conference
#'   on Knowledge Discovery and Data Mining - KDD ’17. New York, New York, USA: ACM Press; 2017. p.
#'   125–34.
#' @examples
#' data <- mp_test_data$train$data[1:1000]
#' window <- 50
#' av <- av_complexity(data, window)
#'
av_complexity <- function(data, window_size, dilution_factor = 0) {
  data <- as.matrix(data)
  data <- as.matrix(data[, 1])
  data_size <- nrow(data)

  if (window_size > data_size / 2) {
    stop("Error: Time series is too short relative to desired window size.", call. = FALSE)
  }
  if (window_size < 4) {
    stop("Error: `window_size` must be at least 4.", call. = FALSE)
  }

  data <- znorm(data)
  profile_size <- data_size - window_size + 1
  av <- matrix(0, profile_size, 1)

  for (j in 1:profile_size) {
    av[j] <- complexity(data[j:(j + window_size - 1)])
  }

  av <- zero_one_norm(av) # zero-one normalize the av

  # Select dilution factor, 0 is no dilution,
  # larger numbers are more dilution
  av <- av + dilution_factor
  av <- av / (dilution_factor + 1)

  av <- list(av = av, w = window_size)
  class(av) <- "AnnotationVector"
  return(av)
}

#' Computes the annotation vector that favors number of zero crossing
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#' @param window_size an `int`. Size of the sliding window.
#'
#' @return Returns the annotation vector for matrix profile correction.
#' @export
#' @references * Dau HA, Keogh E. Matrix Profile V: A Generic Technique to Incorporate Domain
#'   Knowledge into Motif Discovery. In: Proceedings of the 23rd ACM SIGKDD International Conference
#'   on Knowledge Discovery and Data Mining - KDD ’17. New York, New York, USA: ACM Press; 2017. p.
#'   125–34.
#' @examples
#' data <- mp_test_data$train$data[1:1000]
#' window <- 50
#' av <- av_zerocrossing(data, window)
#'
av_zerocrossing <- function(data, window_size) {
  data <- as.matrix(data)
  data <- as.matrix(data[, 1])
  data_size <- nrow(data)

  if (window_size > data_size / 2) {
    stop("Error: Time series is too short relative to desired window size.", call. = FALSE)
  }
  if (window_size < 4) {
    stop("Error: `window_size` must be at least 4.", call. = FALSE)
  }

  data <- znorm(data)
  profile_size <- data_size - window_size + 1
  av <- matrix(0, profile_size, 1)
  for (j in 1:profile_size) {
    av[j] <- zero_crossings(data[j:(j + window_size - 1), ])
  }

  av <- zero_one_norm(av)
  av <- list(av = av, w = window_size)
  class(av) <- "AnnotationVector"
  return(av)
}

#' Computes the annotation vector that suppresses motion artifacts
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#' @param window_size an `int`. Size of the sliding window.
#'
#' @return Returns the annotation vector for matrix profile correction.
#' @export
#' @references * Dau HA, Keogh E. Matrix Profile V: A Generic Technique to Incorporate Domain
#'   Knowledge into Motif Discovery. In: Proceedings of the 23rd ACM SIGKDD International Conference
#'   on Knowledge Discovery and Data Mining - KDD ’17. New York, New York, USA: ACM Press; 2017. p.
#'   125–34.
#' @examples
#' data <- mp_test_data$train$data[1:1000]
#' window <- 50
#' av <- av_motion_artifact(data, window)
#'
av_motion_artifact <- function(data, window_size) {
  data <- as.matrix(data)
  data <- as.matrix(data[, 1])
  data_size <- nrow(data)

  if (window_size > data_size / 2) {
    stop("Error: Time series is too short relative to desired window size.", call. = FALSE)
  }
  if (window_size < 4) {
    stop("Error: `window_size` must be at least 4.", call. = FALSE)
  }

  data <- znorm(data)
  profile_size <- data_size - window_size + 1
  av <- matrix(0, profile_size, 1)

  for (i in 1:profile_size) {
    s <- data[i:(i + window_size - 1), ]
    av[i] <- stats::sd(s)
  }

  cav <- av
  mu <- mean(av)

  cav[av >= mu] <- 0
  cav[av < mu] <- 1

  class(cav) <- "AnnotationVector"
  cav <- list(av = cav, w = window_size)
  return(cav)
}

#' Computes the annotation vector that suppresses stop-word motifs
#'
#' @details
#' The function is intended to be generic. However, its parameters (`stop_word_loc`,
#' `exclusion_zone` and `threshold`) are highly dataset dependent.
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#' @param window_size an `int`. Size of the sliding window.
#' @param stop_word_loc an `int`. The index of stop word location.
#' @param exclusion_zone a `numeric`. Size of the exclusion zone, based on window_size (default is
#'   `1/2`). See details.
#' @param threshold a `numeric`.
#'
#' @return Returns the annotation vector for matrix profile correction.
#' @export
#' @references * Dau HA, Keogh E. Matrix Profile V: A Generic Technique to Incorporate Domain
#'   Knowledge into Motif Discovery. In: Proceedings of the 23rd ACM SIGKDD International Conference
#'   on Knowledge Discovery and Data Mining - KDD ’17. New York, New York, USA: ACM Press; 2017. p.
#'   125–34.
#' @examples
#' data <- mp_test_data$train$data[1:1000]
#' window <- 50
#' av <- av_stop_word(data, window, 150)
#'
av_stop_word <- function(data, window_size, stop_word_loc, exclusion_zone = 1 / 2, threshold = 0.1) {
  data <- as.matrix(data)
  data <- as.matrix(data[, 1])
  data_size <- nrow(data)

  if (window_size > data_size / 2) {
    stop("Error: Time series is too short relative to desired window size.", call. = FALSE)
  }
  if (window_size < 4) {
    stop("Error: `window_size` must be at least 4.", call. = FALSE)
  }
  data <- znorm(data)
  stop_word <- data[stop_word_loc:(stop_word_loc + window_size - 1), ]

  profile_size <- data_size - window_size + 1

  av <- matrix(0, profile_size, 1)

  for (i in 1:profile_size) {
    s <- data[i:(i + window_size - 1), ]
    av[i, ] <- diff2(s, stop_word)
  }

  av <- zero_one_norm(av)

  index <- which(av <= threshold)

  for (i in seq_len(length(index))) {
    if (index[i] < exclusion_zone) {
      av[(index[i] - index[i] + 1):min((index[i] + exclusion_zone - 1), profile_size), ] <- 0
    } else {
      av[(index[i] - exclusion_zone + 1):min((index[i] + exclusion_zone - 1), profile_size), ] <- 0
    }
  }

  av <- list(av = av, w = window_size)
  class(av) <- "AnnotationVector"
  return(av)
}

#' Computes the annotation vector that suppresses hard-limited artifacts
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#' @param window_size an `int`. Size of the sliding window.
#'
#' @return Returns the annotation vector for matrix profile correction.
#' @export
#' @references * Dau HA, Keogh E. Matrix Profile V: A Generic Technique to Incorporate Domain
#'   Knowledge into Motif Discovery. In: Proceedings of the 23rd ACM SIGKDD International Conference
#'   on Knowledge Discovery and Data Mining - KDD ’17. New York, New York, USA: ACM Press; 2017. p.
#'   125–34.
#' @examples
#' data <- mp_test_data$train$data[1:1000]
#' window <- 50
#' av <- av_hardlimit_artifact(data, window)
#'
av_hardlimit_artifact <- function(data, window_size) {
  data <- as.matrix(data)
  data <- as.matrix(data[, 1])
  data_size <- nrow(data)

  if (window_size > data_size / 2) {
    stop("Error: Time series is too short relative to desired window size.", call. = FALSE)
  }
  if (window_size < 4) {
    stop("Error: `window_size` must be at least 4.", call. = FALSE)
  }

  data <- znorm(data)
  max <- max(data)
  min <- min(data)

  profile_size <- data_size - window_size + 1
  av <- matrix(0, profile_size, 1)

  for (i in 1:profile_size) {
    s <- data[i:(i + window_size - 1), ]
    av[i, ] <- length(s[s == max | s == min])
  }

  av <- zero_one_norm(av) # zero-one normalize the av
  av <- 1 - av

  av <- list(av = av, w = window_size)
  class(av) <- "AnnotationVector"
  return(av)
}

#' Corrects the matrix profile using an annotation vector
#'
#' @param .mp The matrix profile.
#' @param annotation_vector The annotation vector.
#'
#' @return Returns the corrected matrix profile
#' @export
#' @references * Dau HA, Keogh E. Matrix Profile V: A Generic Technique to Incorporate Domain
#'   Knowledge into Motif Discovery. In: Proceedings of the 23rd ACM SIGKDD International Conference
#'   on Knowledge Discovery and Data Mining - KDD ’17. New York, New York, USA: ACM Press; 2017. p.
#'   125–34.
#' @examples
#' \dontrun{
#'   av <- av_complexity(data, window)
#'   mpc <- av_apply(mp, av)
#' }
av_apply <- function(.mp, annotation_vector) {
  if (!any(class(.mp) %in% "MatrixProfile")) {
    stop("Error: First argument must be an object of class `MatrixProfile`.")
  }

  if (!(class(annotation_vector) %in% c("AnnotationVector"))) {
    stop("Error: `annotation_vector` must be an object of class `AnnotationVector`.")
  }

  if (.mp$w != annotation_vector$w) {
    warning("Warning: `annotation_vector` window size is not the same as the Matrix Profile.")
  }

  # TODO: do not modify mp
  .mp$mp <- .mp$mp + (1 - annotation_vector$av) * max(.mp$mp)

  return(.mp)
}
