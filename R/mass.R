#' Calculates the distance profile using MASS algorithm
#'
#' Mueen's Algorithm for Similarity Search is The Fastest Similarity Search Algorithm for Time Series Subsequences under Euclidean Distance and Correlation Coefficient.
#'
#' @param data.fft precomputed data product.
#' @param query.window a `vector` of `numeric`. Query window.
#' @param data.size an `int`. The length of the reference data.
#' @param window.size an `int`. Sliding window size.
#' @param data.mean precomputed data moving average.
#' @param data.sd precomputed data moving standard deviation.
#' @param query.mean precomputed query average.
#' @param query.sd precomputed query standard deviation.
#'
#' @return Returns the `distance.profile` for the given query and the `last.product` for STOMP algorithm.
#' @export
#'
#' @seealso [mass.pre()] to precomputation of input values.
#'
#' @references 1. Abdullah Mueen, Yan Zhu, Michael Yeh, Kaveh Kamgar, Krishnamurthy Viswanathan, Chetan Kumar Gupta and Eamonn Keogh (2015), The Fastest Similarity Search Algorithm for Time Series Subsequences under Euclidean Distance
#' @references <https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html>
#'
#' @examples
#' w <- toy_data$sub.len
#' ref.data <- toy_data$data[,1]
#' query.data <- toy_data$data[,1]
#' d.size <- length(ref.data)
#' q.size <- length(query.data)
#'
#' pre <- mass.pre(ref.data, d.size, query.data, q.size, w)
#'
#' dp <- list()
#' for(i in 1:(d.size - w + 1)) {
#'   dp[[i]] <- mass(pre$data.fft, query.data[i:(i-1+w)], d.size, w, pre$data.mean, pre$data.sd,
#'           pre$query.mean[i], pre$query.sd[i])
#' }

mass <- function(data.fft, query.window, data.size, window.size, data.mean, data.sd, query.mean, query.sd) {
  # pre-process query for fft
  query.window <- rev(query.window)
  query.window[(window.size + 1):(window.size + data.size)] <- 0
  # compute the product
  Z <- data.fft * stats::fft(query.window)
  z <- stats::fft(Z, inverse = TRUE) / length(Z)
  # compute the distance profile
  distance.profile <- 2 * (window.size - (z[window.size:data.size] - window.size * data.mean * query.mean) / (data.sd * query.sd))
  last.product <- Re(z[window.size:data.size])

  return(list(distance.profile = distance.profile, last.product = last.product))
}
