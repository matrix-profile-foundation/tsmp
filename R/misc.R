#' Fast implementation of moving standard deviation using filter
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#' @param window.size moving sd window size
#'
#' @return Returns a `vector` with the moving standard deviation
#' @export
#'
#' @examples
#' data.sd <- fast.movsd(toy_data$data[,1], toy_data$sub.len)

fast.movsd <- function(data, window.size) {

  # length of the time series
  data.size <- length(data)

  if (window.size < 2) {
    stop("'window.size' must be at least 2.")
  }

  if (data.size < window.size) {
    stop("'window.size' is too large for this series.")
  }

  # Improve the numerical analysis by subtracting off the series mean
  # this has no effect on the standard deviation.
  data <- data - mean(data)

  # scale the data to have unit variance too. will put that
  # scale factor back into the result at the end
  data.sd <- std(data)
  data <- data / data.sd

  # we will need the squared elements
  data.sqr <- data^2

  b <- matrix(1, 1, window.size)
  s <- sqrt((stats::filter(data.sqr, b, sides = 1) - (stats::filter(data, b, sides = 1)^2) * (1 / window.size)) / (window.size - 1))

  # restore the scale factor that was used before to normalize the data
  s <- s * data.sd
  s <- Re(s)
  s <- s * sqrt((window.size - 1) / window.size)

  return(s[!is.na(s)])
}

#' Fast implementation of moving average using filter
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#' @param window.size moving average window size
#'
#' @return Returns a `vector` with the moving average
#' @export
#' @examples
#' data.avg <- fast.movavg(toy_data$data[,1], toy_data$sub.len)

fast.movavg <- function(data, window.size) {
  data.mean <- stats::filter(data, rep(1 / window.size, window.size), sides = 2)
  return(data.mean[!is.na(data.mean)])
}

#' Population SD, as R always calculate with n-1 (sample), here we fix it
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#'
#' @return Returns the corrected standard deviation from sample to population
#' @keywords internal
#'
std <- function(data) {
  sdx <- stats::sd(data)

  if (sdx == 0) {
    return(sdx)
  }

  return(sqrt((length(data) - 1) / length(data)) * sdx)
}


#' Counts number of zero-crossings
#'
#' Count the number of zero-crossings from the supplied time-domain input vector. A simple method is
#' applied here that can be easily ported to a real-time system that would minimize the number of
#' if-else conditionals.
#'
#' @param data is the input time-domain signal (one dimensional).
#'
#' @return Returns the amount of zero-crossings in the input signal.
#' @author sparafucile17 06/27/04
#' @references <https://www.dsprelated.com/showcode/179.php>
#' @keywords internal
#' @examples
#' count <- zero.crossings(sin(seq(0,10,0.1)))
#'
zero.crossings <- function(data) {
  # initial value
  count <- 0

  data <- as.matrix(data)

  # error checks
  if (length(data) == 1) {
    stop("Input signal must have more than one element")
  }

  if ((ncol(data) != 1) && (nrow(data) != 1)) {
    stop("Input must be one-dimensional")
  }

  # force signal to be a vector oriented in the same direction
  data <- as.vector(data)

  num.samples <- length(data)

  for (i in 2:num.samples) {
    # Any time you multiply to adjacent values that have a sign difference
    # the result will always be negative.  When the signs are identical,
    # the product will always be positive.
    if ((data[i] * data[i - 1]) < 0) {
      count <- count + 1
    }
  }

  return(count)
}

#' Normalizes data for mean Zero and Standard Deviation One
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#'
#' @return Returns the normalized data
#' @keywords internal
#'
znorm <- function(data) {
  data.mean <- mean(data)
  data.dev <- std(data)

  if (is.nan(data.dev) || data.dev <= 0.01) {
    return(data - data.mean)
  }
  else {
    (data - data.mean) / (data.dev)
  }
}

#' Normalizes data between Zero and One
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#'
#' @return Returns the normalized data.
#' @keywords internal
#'
zero.one.norm <- function(data) {
  data <- round(data, 10)

  data <- data - min(data[!is.infinite(data) & !is.na(data)])
  data <- data / max(data[!is.infinite(data) & !is.na(data)])

  return(data)
}

#' Computes the complexity index of the data
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#'
#' @return Returns the complexity index of the data provided (normally a subset)
#' @keywords internal
#' @examples
#' compl <- complexity(sin(seq(0,10,0.1)))
#'
complexity <- function(data) {
  return(sqrt(sum(diff(data)^2)))
}

#' Play sound with `audio`
#'
#' @param data sound data provided by this package
#'
#' @keywords internal
#' @import audio
#'
beep <- function(data) {
  if (!(is.null(audio::audio.drivers()) || nrow(audio::audio.drivers()) == 0)) {
    tryCatch({
      audio::play(data)
    },
    error = function(cond) {
      message("Failed to play audio alert")
      message(cond)
    },
    warning = function(cond) {
      message("Something went wrong playing audio alert")
      message(cond)
    }
    )
  }
  Sys.sleep(1)
  invisible()
}

#' Distance between two matrices
#'
#' Computes the Euclidean distance between rows of two matrices.
#'
#' @param X a `matrix`.
#' @param Y a `matrix`.
#'
#' @return Returns a `matrix` of size m x n if x is of size m x k and y is of size n x k.
#' @keywords internal

diff2 <- function(X, Y) {
  if (!is.numeric(X) || !is.numeric(Y)) {
    stop("X and Y must be numeric vectors or matrices.")
  }
  if (is.vector(X)) {
    dim(X) <- c(1, length(X))
  }
  if (is.vector(Y)) {
    dim(Y) <- c(1, length(Y))
  }
  if (ncol(X) != ncol(Y)) {
    stop("X and Y must have the same number of columns.")
  }
  m <- nrow(X)
  n <- nrow(Y)
  XY <- X %*% t(Y)
  XX <- matrix(rep(apply(X * X, 1, sum), n), m, n, byrow = F)
  YY <- matrix(rep(apply(Y * Y, 1, sum), m), m, n, byrow = T)
  sqrt(pmax(XX + YY - 2 * XY, 0))
}

#' Global constants
#'
#' @return Returns a `list` with the global constants
#' @keywords internal
#'
vars <- function() {
  eps <- .Machine$double.eps^0.5

  return(list(eps = eps))
}
