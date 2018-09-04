#' Calculates the distance profile using MASS algorithm
#'
#' Mueen's Algorithm for Similarity Search is The Fastest Similarity Search Algorithm for Time
#' Series Subsequences under Euclidean Distance and Correlation Coefficient.
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
#' @return Returns the `distance.profile` for the given query and the `last.product` for STOMP
#'   algorithm.
#' @export
#'
#' @seealso [mass.pre()] to precomputation of input values.
#'
#' @references * Abdullah Mueen, Yan Zhu, Michael Yeh, Kaveh Kamgar, Krishnamurthy Viswanathan,
#'   Chetan Kumar Gupta and Eamonn Keogh (2015), The Fastest Similarity Search Algorithm for Time
#'   Series Subsequences under Euclidean Distance
#' @references Website: <https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html>
#'
#' @examples
#' w <- mp_toy_data$sub.len
#' ref.data <- mp_toy_data$data[,1]
#' query.data <- mp_toy_data$data[,1]
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

mass <- function(data_fft, query_window, data_size, window_size, data_mean, data_sd, query_mean, query_sd) {
  # pre-process query for fft
  query_window <- rev(query_window)
  query_window[(window_size + 1):(window_size + data_size)] <- 0
  # compute the product
  prod <- data_fft * stats::fft(query_window)
  z <- stats::fft(prod, inverse = TRUE) / length(prod)
  # compute the distance profile
  distance_profile <- 2 * (window_size - (z[window_size:data_size] - window_size * data_mean * query_mean) / (data_sd * query_sd))
  last_product <- z[window_size:data_size]

  return(list(distance_profile = distance_profile, last_product = last_product))
}
