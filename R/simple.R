#' Compute the similarity join for Sound data.
#'
#' Compute the similarity join for Sound data.
#'
#' `verbose` changes how much information is printed by this function; `0` means nothing, `1` means text, `2` means text and sound.
#'
#' @param data a `matrix` of `numeric`, where each colums is a time series. Accepts `list` and `data.frame` too.
#' @param window.size an `int` with the size of the sliding window.
#' @param exclusion.zone a `numeric`. Size of the exclusion zone, based on query size (default is `1/2`).
#' @param verbose an `int`. See details. (Default is `2`).
#'
#' @return Returns a list with the Matrix Profile `mp` and Profile Index `pi`.
#'
#' @export
#' @references 1. Silva D, Yeh C, Batista G, Keogh E. Simple: Assessing Music Similarity Using Subsequences Joins. Proc 17th ISMIR Conf. 2016;23–30.
#' @references 2. Silva DF, Yeh C-CM, Zhu Y, Batista G, Keogh E. Fast Similarity Matrix Profile for Music Analysis and Exploration. IEEE Trans Multimed. 2018;14(8):1–1.
#' @references Website: <https://sites.google.com/view/simple-fast>
#' @references Website: <https://sites.google.com/site/ismir2016simple/home>
#'
#' @examples
#' w <- 30
#' data <- toy_data$data # 3 dimensions matrix
#' result <- simple.fast(data, w, verbose = 0)
#'
simple.fast <- function(data, window.size, exclusion.zone = 1 / 2, verbose = 2) {
  ## get various length
  exclusion.zone <- floor(window.size * exclusion.zone)

  ## transform data list into matrix
  if (is.list(data)) {
    data.size <- length(data[[1]])
    n.dim <- length(data)

    for (i in 1:n.dim) {
      len <- length(data[[i]])
      # Fix TS size with NaN
      if (len < data.size) {
        data[[i]] <- c(data[[i]], rep(NA, data.size - len))
      }
    }
    # transform data into matrix (each column is a TS)
    data <- sapply(data, cbind)
  } else if (is.matrix(data) || is.data.frame(data)) {
    if (is.data.frame(data)) {
      data <- as.matrix(data)
    } # just to be uniform
    if (ncol(data) > nrow(data)) {
      data <- t(data)
    }
    data.size <- nrow(data)
    n.dim <- ncol(data)
  } else if (is.vector(data)) {
    data.size <- length(data)
    n.dim <- 1
    # transform data into 1-col matrix
    data <- as.matrix(data) # just to be uniform
  } else {
    stop("Unknown type of data. Must be: matrix, data.frame, vector or list")
  }

  ## check input
  if (window.size > data.size / 2) {
    stop("Error: Time series is too short relative to desired subsequence length")
  }
  if (window.size < 4) {
    stop("Error: Subsequence length must be at least 4")
  }

  ## initialization
  matrix.profile.size <- data.size - window.size + 1
  matrix.profile <- rep(Inf, matrix.profile.size)
  profile.index <- rep(0, matrix.profile.size)

  if (verbose > 0) {
    pb <- utils::txtProgressBar(min = 0, max = matrix.profile.size, style = 3, width = 80)
    on.exit(close(pb))
  }
  if (verbose > 1) {
    on.exit(audio::play(sounds[[1]]), TRUE)
  }

  ## compute necessary values
  res <- mass.simple.pre(data, data.size, window.size = window.size)
  data.fft <- res$data.fft
  sumx2 <- res$sumx2

  ## compute first distance profile
  query.window <- data[1:window.size, ]
  res <- mass.simple(data.fft, query.window, data.size, window.size, sumx2)
  distance.profile <- res$distance.profile
  first.product <- last.product <- res$last.product
  sumy2 <- res$sumy2
  dropval <- query.window[1, ]
  distance.profile[1:exclusion.zone] <- Inf


  ind <- which.min(distance.profile)
  profile.index[1] <- ind
  matrix.profile[1] <- distance.profile[ind]

  ## compute the remainder of the matrix profile
  for (i in 2:matrix.profile.size) {

    # compute the distance profile
    if (verbose > 0) {
      utils::setTxtProgressBar(pb, i)
    }

    query.window <- data[i:(i + window.size - 1), ]

    sumy2 <- sumy2 - dropval^2 + query.window[window.size, ]^2

    for (j in 1:n.dim) {
      last.product[2:(data.size - window.size + 1), j] <- last.product[1:(data.size - window.size), j] -
        data[1:(data.size - window.size), j] * dropval[j] +
        data[(window.size + 1):data.size, j] * query.window[window.size, j]
    }

    last.product[1, ] <- first.product[i, ]
    dropval <- query.window[1, ]

    distance.profile <- matrix(0, nrow(sumx2), 1)

    for (j in 1:n.dim) {
      distance.profile <- distance.profile + sumx2[, j] - 2 * last.product[, j] + sumy2[j]
    }

    exc.st <- max(1, i - exclusion.zone)
    exc.ed <- min(matrix.profile.size, i + exclusion.zone)
    distance.profile[exc.st:exc.ed] <- Inf

    ind <- which.min(distance.profile)
    profile.index[i] <- ind
    matrix.profile[i] <- distance.profile[ind]
  }

  return(list(mp = matrix.profile, pi = profile.index))
}

#' Precomputes several values used on MASS
#'
#' The difference of this function to [mass.pre()] is that this does not normalize data. Specific for this domain.
#'
#' @param data a `matrix` of `numeric`. Reference Time Series.
#' @param data.size an `int`. Reference Time Series size.
#' @param window.size an `int`. Sliding window size.
#'
#' @return Returns `data.fft` and `sumx2`.
#' @keywords internal
#'
#' @references Abdullah Mueen, Yan Zhu, Michael Yeh, Kaveh Kamgar, Krishnamurthy Viswanathan, Chetan Kumar Gupta and Eamonn Keogh (2015), The Fastest Similarity Search Algorithm for Time Series Subsequences under Euclidean Distance.
#' @references <https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html>


mass.simple.pre <- function(data, data.size, window.size) {
  if (nrow(data) < ncol(data)) {
    data <- t(data)
  }

  n.dim <- ncol(data)

  data <- rbind(data, matrix(0, data.size, n.dim))

  data.fft <- apply(data, 2, stats::fft)
  cum_sumx2 <- apply(data^2, 2, cumsum)

  sumx2 <- cum_sumx2[window.size:data.size, ] - rbind(rep(0, n.dim), cum_sumx2[1:(data.size - window.size), ])

  return(list(data.fft = data.fft, sumx2 = sumx2))
}

#' Calculates the distance profile using MASS algorithm
#'
#' Mueen's Algorithm for Similarity Search is The Fastest Similarity Search Algorithm for Time Series Subsequences under Euclidean Distance and Correlation Coefficient.
#' The difference of this function to [mass()] is that this does not normalize data. Specific for this domain.
#'
#' @param data.fft precomputed data product.
#' @param query.window a `matrix` of `numeric`. Query window.
#' @param data.size an `int`. The length of the reference data.
#' @param window.size an `int`. Sliding window size.
#' @param sumx2 precomputed sum of squares
#'
#' @return Returns the `distance.profile` for the given query and the `last.product` for STOMP algorithm and `sumy2`.
#' @keywords internal
#'
#' @references Abdullah Mueen, Yan Zhu, Michael Yeh, Kaveh Kamgar, Krishnamurthy Viswanathan, Chetan Kumar Gupta and Eamonn Keogh (2015), The Fastest Similarity Search Algorithm for Time Series Subsequences under Euclidean Distance
#' @references <https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html>
#'


mass.simple <- function(data.fft, query.window, data.size, window.size, sumx2) {
  if (nrow(data.fft) < ncol(data.fft)) {
    data.fft <- t(data.fft)
  }

  n.dim <- ncol(data.fft)

  if (ncol(query.window) != n.dim) {
    query.window <- t(query.window)
  }

  # pre-process query for fft
  query.window <- apply(query.window, 2, rev)
  query.window <- rbind(query.window, matrix(0, 2 * data.size - window.size, n.dim))

  query.fft <- apply(query.window, 2, stats::fft)
  # compute the product
  Z <- data.fft * query.fft
  z <- apply(Z, 2, function(x){
    stats::fft(x, inverse = TRUE) / length(x)
  })

  sumy2 <- apply(query.window^2, 2, sum)

  last.product <- Re(z[window.size:data.size, ])

  distance.profile <- matrix(0, nrow(sumx2), 1)

  for (i in 1:n.dim) {
    distance.profile <- distance.profile + sumx2[, i] - 2 * last.product[, i] + sumy2[i]
  }

  return(list(distance.profile = distance.profile, last.product = last.product, sumy2 = sumy2))
}
