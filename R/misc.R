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

fast_movsd <- function(data, window_size) {

  # length of the time series
  data_size <- length(data)

  if (window_size < 2) {
    stop("Error: 'window_size' must be at least 2.", call. = FALSE)
  }

  if (data_size < window_size) {
    stop("Error: 'window_size' is too large for this series.", call. = FALSE)
  }

  # Improve the numerical analysis by subtracting off the series mean
  # this has no effect on the standard deviation.
  data <- data - mean(data)

  # scale the data to have unit variance too. will put that
  # scale factor back into the result at the end
  data_sd <- std(data)
  data <- data / data_sd

  # we will need the squared elements
  data_sqr <- data^2

  b <- matrix(1, 1, window_size)
  s <- sqrt((stats::filter(data_sqr, b, sides = 1) - (stats::filter(data, b, sides = 1)^2) * (1 / window_size)) / (window_size - 1))

  # restore the scale factor that was used before to normalize the data
  s <- s * data_sd
  s <- Re(s)
  s <- s * sqrt((window_size - 1) / window_size)

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

fast_movavg <- function(data, window_size) {
  data_mean <- stats::filter(data, rep(1 / window_size, window_size), sides = 2)
  return(data_mean[!is.na(data_mean)])
}

#' Population SD, as R always calculate with n-1 (sample), here we fix it
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#'
#' @return Returns the corrected standard deviation from sample to population
#' @keywords internal
#' @noRd
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
#' @noRd
#'
zero_crossings <- function(data) {
  # initial value
  count <- 0

  data <- as.matrix(data)

  # error checks
  if (length(data) == 1) {
    stop("Error: Input signal must have more than one element.", call. = FALSE)
  }

  if ((ncol(data) != 1) && (nrow(data) != 1)) {
    stop("Error: Input must be one-dimensional.", call. = FALSE)
  }

  # force signal to be a vector oriented in the same direction
  data <- as.vector(data)

  num_samples <- length(data)

  for (i in 2:num_samples) {
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
#' @noRd
#'
znorm <- function(data) {
  data_mean <- mean(data)
  data_dev <- std(data)

  if (is.nan(data_dev) || data_dev <= 0.01) {
    return(data - data_mean)
  }
  else {
    (data - data_mean) / (data_dev)
  }
}

#' Normalizes data between Zero and One
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#'
#' @return Returns the normalized data.
#' @keywords internal
#' @noRd
#'
zero_one_norm <- function(data) {
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
#' @noRd
#'
complexity <- function(data) {
  return(sqrt(sum(diff(data)^2)))
}

#' Play sound with `audio`
#'
#' @param data sound data provided by this package
#'
#' @keywords internal
#' @noRd
#' @import audio
#'
beep <- function(data) {
  if (!(is.null(audio::audio.drivers()) || nrow(audio::audio.drivers()) == 0)) {
    tryCatch({
      audio::play(data)
    },
    error = function(cond) {
      message("Warning: Failed to play audio alert")
      message(cond)
    },
    warning = function(cond) {
      message("Warning: Something went wrong playing audio alert")
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
#' @param x a `matrix`.
#' @param y a `matrix`.
#'
#' @return Returns a `matrix` of size m x n if x is of size m x k and y is of size n x k.
#' @keywords internal
#' @noRd

diff2 <- function(x, y) {
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("Error: `x` and `y` must be numeric vectors or matrices.", call. = FALSE)
  }
  if (is.vector(x)) {
    dim(x) <- c(1, length(x))
  }
  if (is.vector(y)) {
    dim(y) <- c(1, length(y))
  }
  if (ncol(x) != ncol(y)) {
    stop("Error: `x` and `y` must have the same number of columns.", call. = FALSE)
  }
  m <- nrow(x)
  n <- nrow(y)
  xy <- x %*% t(y)
  xx <- matrix(rep(apply(x * x, 1, sum), n), m, n, byrow = FALSE)
  yy <- matrix(rep(apply(y * y, 1, sum), m), m, n, byrow = TRUE)
  sqrt(pmax(xx + yy - 2 * xy, 0))
}

#' Global constants
#'
#' @return Returns a `list` with the global constants
#' @keywords internal
#' @noRd
#'
vars <- function() {
  eps <- .Machine$double.eps^0.5

  return(list(eps = eps))
}
