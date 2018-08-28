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
#' av <- av.complexity(data, window)
#'
av.complexity <- function(data, window.size, dilution.factor = 0) {
  data <- znorm(data) # data is a row vector
  profile.size <- length(data) - window.size + 1
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
#' av <- av.zerocrossing(data, window)
#'
av.zerocrossing <- function(data, window.size) {
  # data is a row vector
  data <- znorm(data)
  profile.size <- length(data) - window.size + 1
  av <- matrix(0, profile.size, 1)
  for (j in 1:profile.size) {
    av[j] <- zero.crossings(data[j:(j + window.size - 1)])
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
#' av <- av.motion.artifact(data, window)
#'
av.motion.artifact <- function(data, window.size) {
  data <- znorm(data)
  profile.size <- length(data) - window.size + 1
  av <- matrix(0, profile.size, 1)

  for (i in 1:profile.size) {
    s <- data[i:(i + window.size - 1)]
    av[i] <- sd(s)
  }

  cav <- av
  mu <- mean(av)

  cav[av >= mu] <- 0
  cav[av < mu] <- 1

  return(cav)
}

# The function is intended to be generic. However, its parameters
# (threshold, exclusion.zone and stop.word.location) are dataset-dependent
# The parameters' default value are for specifically for dataset
# ECG_LTAF-71.mat. Recommened window.size of this dataset is 150.
#
##
#' Computes the annotation vector that supresses stop-word motifs
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
#' av <- av.stop.word(data, window)
#'
av.stop.word <- function(data, window.size) {
  ## TODO: NOT WORKING
  # the following parameters are dataset-dependent
  threshold <- 0.1
  exclusion.zone <- 450
  stop.word.location <- 63

  data <- znorm(data)
  stop.word <- data[stop.word.location:(stop.word.location + window.size - 1)]

  profile.size <- length(data) - window.size + 1

  av <- matrix(0, profile.size, 1)

  for (i in 1:profile.size) {
    s <- data[i:(i + window.size - 1)]
    av[i] <- pdist2(s, stop.word)
  }

  av <- zero.one.norm(av)

  index <- which(av <= threshold)

  for (i in 1:length(index)) {
    if (index[i] < exclusion.zone) {
      av[(index[i] - index[i] + 1):(index[i] + exclusion.zone - 1)] <- 0
    } else {
      av[(index[i] - exclusion.zone + 1):(index[i] + exclusion.zone - 1)] <- 0
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
#' av <- av.hardlimit.artifact(data, window)
#'
av.hardlimit.artifact <- function(data, window.size) {
  data <- znorm(data)
  max <- max(data)
  min <- min(data)

  profile.size <- length(data) - window.size + 1
  av <- matrix(0, profile.size, 1)

  for (i in 1:profile.size) {
    s <- data[i:(i + window.size - 1)]
    av[i] <- length(s[s == max | s == min])
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
#' av <- av.complexity(data, window)
#' mpc <- av.apply(mp, av)
#'
av.apply <- function(matrix.profile, annotation.vector) {
  corrected.mp <- matrix.profile + (1 - annotation.vector) * max(matrix.profile)

  return(corrected.mp)
}

# Guided motif search
