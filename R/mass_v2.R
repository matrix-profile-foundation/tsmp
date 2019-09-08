#' Calculates the distance profile using MASS_V2 algorithm
#'
#' Mueen's Algorithm for Similarity Search is The Fastest Similarity Search Algorithm for Time
#' Series Subsequences under Euclidean Distance and Correlation Coefficient.
#'
#' @param data_fft precomputed data product.
#' @param query_window a `vector` of `numeric`. Query window.
#' @param data_size an `int`. The length of the reference data.
#' @param window_size an `int`. Sliding window size.
#' @param data_mean precomputed data moving average.
#' @param data_sd precomputed data moving standard deviation.
#' @param query_mean precomputed query average.
#' @param query_sd precomputed query standard deviation.
#' @param \dots just a placeholder to catch unused parameters.
#'
#' @return Returns the `distance_profile` for the given query and the `last_product` for STOMP
#'   algorithm.
#'
#' @seealso [mass_pre()] to precomputation of input values.
#'
#' @references * Abdullah Mueen, Yan Zhu, Michael Yeh, Kaveh Kamgar, Krishnamurthy Viswanathan,
#'   Chetan Kumar Gupta and Eamonn Keogh (2015), The Fastest Similarity Search Algorithm for Time
#'   Series Subsequences under Euclidean Distance
#' @references Website: <https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html>
#'
#' @keywords internal
#'
#' @examples
#' w <- mp_toy_data$sub_len
#' ref_data <- mp_toy_data$data[, 1]
#' query_data <- mp_toy_data$data[, 1]
#' d_size <- length(ref_data)
#' q_size <- length(query_data)
#'
#' pre <- tsmp:::mass_pre(ref_data, query_data, w)
#'
#' dp <- list()
#' for (i in 1:(d_size - w + 1)) {
#'   dp[[i]] <- tsmp:::mass_v2(
#'     query_data[i:(i - 1 + w)], pre$window_size,
#'     pre$data_fft, pre$data_size, pre$data_mean, pre$data_sd,
#'     pre$query_mean[i], pre$query_sd[i]
#'   )
#' }
mass_v2 <- function(query_window, window_size, data_fft, data_size, data_mean, data_sd, query_mean, query_sd, ...) {
  # pre-process query for fft
  query_window <- rev(query_window)
  pad_size <- length(data_fft)
  query_window[(window_size + 1):pad_size] <- 0
  # compute the product
  prod <- data_fft * stats::fft(query_window)
  z <- Re(stats::fft(prod, inverse = TRUE) / length(prod))
  # compute the distance profile
  last_product <- z[window_size:data_size]
  distance_profile <- 2 * (window_size - (last_product - window_size * data_mean * query_mean) / (data_sd * query_sd))
  distance_profile[distance_profile < 0] <- 0

  return(list(distance_profile = distance_profile, last_product = last_product))
}
