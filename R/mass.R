#' Calculates the distance profile using MASS algorithm
#'
#' Mueen's Algorithm for Similarity Search is The Fastest Similarity Search Algorithm for Time Series Subsequences under Euclidean Distance and Correlation Coefficient
#'
#' @param data.fft precomputed data product
#' @param query.window query window vector
#' @param data.size length of data
#' @param window.size length of query moving window
#' @param data.mean precomputed data moving average
#' @param data.sd precomputed data moving standard deviation
#' @param query.mean precomputed query average
#' @param query.sd precomputed query standard deviation
#'
#' @return Returns the distance profile for the given query and the last computed product
#'
#' @references Abdullah Mueen, Yan Zhu, Michael Yeh, Kaveh Kamgar, Krishnamurthy Viswanathan, Chetan Kumar Gupta and Eamonn Keogh (2015), The Fastest Similarity Search Algorithm for Time Series Subsequences under Euclidean Distance
#' @references <https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html>

mass <- function(data.fft, query.window, data.size, window.size, data.mean, data.sd, query.mean, query.sd) {
  # pre-process query for fft
  query.window <- rev(query.window)
  query.window[(window.size + 1):(window.size + data.size)] <- 0
  # compute the product
  Z <- data.fft * fft(query.window)
  z <- fft(Z, inverse = TRUE) / length(Z)
  # compute the distance profile
  distance.profile <- 2 * (window.size - (z[window.size:data.size] - window.size * data.mean * query.mean) / (data.sd * query.sd))
  last.product <- Re(z[window.size:data.size])

  return(list(distance.profile = distance.profile, last.product = last.product))
}
