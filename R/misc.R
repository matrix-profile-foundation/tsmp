#' Fast implementation of moving standard deviation using filter
#'
#' @param data a vector
#' @param n moving sd window
#'
#' @return Returns a vector with the moving standard deviation
#' @export

fast.movsd <- function(data, n) {

  # length of the time series
  data.size <- length(data)

  if (n < 2) {
    stop("'n' must be at least 2.")
  }

  if (data.size < n) {
    stop("'n' is too large for this series.")
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

  B <- matrix(1, 1, n)
  s <- sqrt((stats::filter(data.sqr, B, sides = 1) - (stats::filter(data, B, sides = 1)^2) * (1 / n)) / (n - 1))

  # restore the scale factor that was used before to normalize the data
  s <- s * data.sd
  s <- Re(s)
  s <- s * sqrt((n - 1) / n)

  return(s[!is.na(s)])
}

#' Fast implementation of moving average using filter
#'
#' @param data a vector
#' @param n moving average window
#'
#' @return Returns a vector with the moving average
#' @export
#'

fast.movavg <- function(data, n) {
  data.mean <- stats::filter(data, rep(1 / n, n), sides = 2)
  return(data.mean[!is.na(data.mean)])
}

#' Population SD, as R always calculate with n-1 (sample), here we fix it
#'
#' @param x a vector
#'
#' @return Returns a corrected standard deviation from sample to population
#' @keywords internal
#'
std <- function(x) {
  sdx <- sd(x)

  if (sdx == 0)
    return(sdx)

  return(sqrt((length(x) - 1) / length(x)) * sdx)
}
