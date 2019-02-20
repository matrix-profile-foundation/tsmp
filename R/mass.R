#' Calculates the distance profile using MASS algorithm
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
#'
#' @return Returns the `distance_profile` for the given query and the `last_product` for STOMP
#'   algorithm.
#' @export
#'
#' @seealso [mass_pre()] to precomputation of input values.
#'
#' @references * Abdullah Mueen, Yan Zhu, Michael Yeh, Kaveh Kamgar, Krishnamurthy Viswanathan,
#'   Chetan Kumar Gupta and Eamonn Keogh (2015), The Fastest Similarity Search Algorithm for Time
#'   Series Subsequences under Euclidean Distance
#' @references Website: <https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html>
#'
#' @examples
#' w <- mp_toy_data$sub_len
#' ref_data <- mp_toy_data$data[, 1]
#' query_data <- mp_toy_data$data[, 1]
#' d_size <- length(ref_data)
#' q_size <- length(query_data)
#' 
#' pre <- mass_pre(ref_data, d_size, query_data, q_size, w)
#' 
#' dp <- list()
#' for (i in 1:(d_size - w + 1)) {
#'   dp[[i]] <- mass(
#'     pre$data_fft, query_data[i:(i - 1 + w)], d_size, w, pre$data_mean, pre$data_sd,
#'     pre$query_mean[i], pre$query_sd[i]
#'   )
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
