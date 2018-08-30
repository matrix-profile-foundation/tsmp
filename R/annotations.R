#' Computes the annotation vector that favors complexity
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#' @param window.size an `int`. Size of the sliding window.
#' @param dilution.factor a `numeric`. (Default is `0`). Larger numbers means more dilution.
#'
#' @return Returns the annotation vector for matrix profile correction.
#' @export
#' @references * Dau HA, Keogh E. Matrix Profile V: A Generic Technique to Incorporate Domain
#'   Knowledge into Motif Discovery. In: Proceedings of the 23rd ACM SIGKDD International Conference
#'   on Knowledge Discovery and Data Mining - KDD ’17. New York, New York, USA: ACM Press; 2017. p.
#'   125–34.
#' @examples
#' data <- test_data$train$data[1:1000]
#' window <- 50
#' av <- av.complexity(data, window)
#'
av.complexity <- function(data, window.size, dilution.factor = 0) {
  data <- as.matrix(data)
  data <- as.matrix(data[, 1])
  data.size <- nrow(data)

  if (window.size > data.size / 2) {
    stop("Error: Time series is too short relative to desired window size")
  }
  if (window.size < 4) {
    stop("Error: Window size must be at least 4")
  }

  data <- znorm(data)
  profile.size <- data.size - window.size + 1
  av <- matrix(0, profile.size, 1)

  for (j in 1:profile.size) {
    av[j] <- complexity(data[j:(j + window.size - 1)])
  }

  av <- zero.one.norm(av) # zero-one normalize the av

  # Select dilution factor, 0 is no dilution,
  # larger numbers are more dilution
  av <- av + dilution.factor
  av <- av / (dilution.factor + 1)

  return(av)
}

#' Computes the annotation vector that favors number of zero crossing
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#' @param window.size an `int`. Size of the sliding window.
#'
#' @return Returns the annotation vector for matrix profile correction.
#' @export
#' @references * Dau HA, Keogh E. Matrix Profile V: A Generic Technique to Incorporate Domain
#'   Knowledge into Motif Discovery. In: Proceedings of the 23rd ACM SIGKDD International Conference
#'   on Knowledge Discovery and Data Mining - KDD ’17. New York, New York, USA: ACM Press; 2017. p.
#'   125–34.
#' @examples
#' data <- test_data$train$data[1:1000]
#' window <- 50
#' av <- av.zerocrossing(data, window)
#'
av.zerocrossing <- function(data, window.size) {
  data <- as.matrix(data)
  data <- as.matrix(data[, 1])
  data.size <- nrow(data)

  if (window.size > data.size / 2) {
    stop("Error: Time series is too short relative to desired window size")
  }
  if (window.size < 4) {
    stop("Error: Window size must be at least 4")
  }

  data <- znorm(data)
  profile.size <- data.size - window.size + 1
  av <- matrix(0, profile.size, 1)
  for (j in 1:profile.size) {
    av[j] <- zero.crossings(data[j:(j + window.size - 1), ])
  }

  av <- zero.one.norm(av)

  return(av)
}

#' Computes the annotation vector that suppresses motion artifacts
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#' @param window.size an `int`. Size of the sliding window.
#'
#' @return Returns the annotation vector for matrix profile correction.
#' @export
#' @references * Dau HA, Keogh E. Matrix Profile V: A Generic Technique to Incorporate Domain
#'   Knowledge into Motif Discovery. In: Proceedings of the 23rd ACM SIGKDD International Conference
#'   on Knowledge Discovery and Data Mining - KDD ’17. New York, New York, USA: ACM Press; 2017. p.
#'   125–34.
#' @examples
#' data <- test_data$train$data[1:1000]
#' window <- 50
#' av <- av.motion.artifact(data, window)
#'
av.motion.artifact <- function(data, window.size) {
  data <- as.matrix(data)
  data <- as.matrix(data[, 1])
  data.size <- nrow(data)

  if (window.size > data.size / 2) {
    stop("Error: Time series is too short relative to desired window size")
  }
  if (window.size < 4) {
    stop("Error: Window size must be at least 4")
  }

  data <- znorm(data)
  profile.size <- data.size - window.size + 1
  av <- matrix(0, profile.size, 1)

  for (i in 1:profile.size) {
    s <- data[i:(i + window.size - 1), ]
    av[i] <- stats::sd(s)
  }

  cav <- av
  mu <- mean(av)

  cav[av >= mu] <- 0
  cav[av < mu] <- 1

  return(cav)
}

#' Computes the annotation vector that suppresses stop-word motifs.
#'
#' Computes the annotation vector that suppresses stop-word motifs.
#'
#' The function is intended to be generic. However, its parameters (`stop.word.loc`,
#' `exclusion.zone` and `threshold`) are highly dataset dependant.
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#' @param window.size an `int`. Size of the sliding window.
#' @param stop.word.loc an `int`. The index of stop word location.
#' @param exclusion.zone a `numeric`. Size of the exclusion zone, based on window.size (default is
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
#' data <- test_data$train$data[1:1000]
#' window <- 50
#' av <- av.stop.word(data, window, 150)
#'
av.stop.word <- function(data, window.size, stop.word.loc, exclusion.zone = 1 / 2, threshold = 0.1) {
  data <- as.matrix(data)
  data <- as.matrix(data[, 1])
  data.size <- nrow(data)

  if (window.size > data.size / 2) {
    stop("Error: Time series is too short relative to desired window size")
  }
  if (window.size < 4) {
    stop("Error: Window size must be at least 4")
  }
  data <- znorm(data)
  stop.word <- data[stop.word.loc:(stop.word.loc + window.size - 1), ]

  profile.size <- data.size - window.size + 1

  av <- matrix(0, profile.size, 1)

  for (i in 1:profile.size) {
    s <- data[i:(i + window.size - 1), ]
    av[i, ] <- diff2(s, stop.word)
  }

  av <- zero.one.norm(av)

  index <- which(av <= threshold)

  for (i in 1:length(index)) {
    if (index[i] < exclusion.zone) {
      av[(index[i] - index[i] + 1):min((index[i] + exclusion.zone - 1), profile.size), ] <- 0
    } else {
      av[(index[i] - exclusion.zone + 1):min((index[i] + exclusion.zone - 1), profile.size), ] <- 0
    }
  }

  return(av)
}

#' Computes the annotation vector that suppresses hard-limited artifacts
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#' @param window.size an `int`. Size of the sliding window.
#'
#' @return Returns the annotation vector for matrix profile correction.
#' @export
#' @references * Dau HA, Keogh E. Matrix Profile V: A Generic Technique to Incorporate Domain
#'   Knowledge into Motif Discovery. In: Proceedings of the 23rd ACM SIGKDD International Conference
#'   on Knowledge Discovery and Data Mining - KDD ’17. New York, New York, USA: ACM Press; 2017. p.
#'   125–34.
#' @examples
#' data <- test_data$train$data[1:1000]
#' window <- 50
#' av <- av.hardlimit.artifact(data, window)
#'
av.hardlimit.artifact <- function(data, window.size) {
  data <- as.matrix(data)
  data <- as.matrix(data[, 1])
  data.size <- nrow(data)

  if (window.size > data.size / 2) {
    stop("Error: Time series is too short relative to desired window size")
  }
  if (window.size < 4) {
    stop("Error: Window size must be at least 4")
  }

  data <- znorm(data)
  max <- max(data)
  min <- min(data)

  profile.size <- data.size - window.size + 1
  av <- matrix(0, profile.size, 1)

  for (i in 1:profile.size) {
    s <- data[i:(i + window.size - 1), ]
    av[i, ] <- length(s[s == max | s == min])
  }

  av <- zero.one.norm(av) # zero-one normalize the av
  av <- 1 - av

  return(av)
}

#' Corrects the matrix profile using an annotation vector
#'
#' @param matrix.profile The matrix profile.
#' @param annotation.vector The annotation vector.
#'
#' @return Returns the corrected matrix profile
#' @export
#' @references * Dau HA, Keogh E. Matrix Profile V: A Generic Technique to Incorporate Domain
#'   Knowledge into Motif Discovery. In: Proceedings of the 23rd ACM SIGKDD International Conference
#'   on Knowledge Discovery and Data Mining - KDD ’17. New York, New York, USA: ACM Press; 2017. p.
#'   125–34.
#' @examples
#' \dontrun{
#'   av <- av.complexity(data, window)
#'   mpc <- av.apply(mp, av)
#' }
av.apply <- function(matrix.profile, annotation.vector) {
  corrected.mp <- matrix.profile + (1 - annotation.vector) * max(matrix.profile)

  return(corrected.mp)
}

# Guided motif search
