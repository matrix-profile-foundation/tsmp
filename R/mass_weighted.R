#' Calculates the distance profile using MASS_WEIGHTED algorithm
#'
#' Mueen's Algorithm for Similarity Search is The Fastest Similarity Search Algorithm for Time
#' Series Subsequences under Euclidean Distance and Correlation Coefficient.
#'
#' @param data_fft precomputed data product.
#' @param query_window a `vector` of `numeric`. Query window.
#' @param data_size an `int`. The length of the reference data.
#' @param window_size an `int`. Sliding window size.
#' @param data_mean precomputed data moving average.
#' @param data_pre precomputed weighted data product.
#' @param weight a `vector` of `numeric` with the same length of the `window_size`.
#' @param data_sd precomputed data moving standard deviation.
#' @param \dots just a placeholder to catch unused parameters.
#'
#' @return Returns the `distance_profile` for the given query and the `last_product` for STOMP
#'   algorithm.
#' @keywords internal
#'
#' @seealso [mass_pre_w()] to precomputation of input values.
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
#' weight <- c(
#'   rep(1, mp_toy_data$sub_len / 3), rep(0.5, mp_toy_data$sub_len / 3),
#'   rep(1, mp_toy_data$sub_len / 3)
#' )
#'
#' pre <- tsmp:::mass_pre_w(ref_data, query_data, w, weight)
#'
#' dp <- list()
#' for (i in 1:(pre$data_size - w + 1)) {
#'   dp[[i]] <- tsmp:::mass_weighted(
#'     query_data[i:(i - 1 + w)], pre$window_size, pre$data_fft, pre$data_size,
#'     pre$data_mean, pre$data_sd, pre$data_pre, pre$weight
#'   )
#' }
mass_weighted <- function(query_window, window_size, data_fft, data_size, data_mean, data_sd,
                          data_pre, weight, ...) {
  if (length(weight) != window_size) {
    stop("'weight' must have the same length as the 'window_size'.")
  }

  # normalized query
  query_window <- (query_window - mean(query_window)) / std(query_window)

  sumwy <- sum(weight * query_window)
  sumwy2 <- sum(weight * (query_window^2))

  query_window <- rev(query_window)
  query_window[(window_size + 1):data_size] <- 0
  weight <- rev(weight)
  weight[(window_size + 1):data_size] <- 0

  queryw_fft <- stats::fft(weight * query_window)
  data_fft_queryw_fft <- data_fft * queryw_fft

  data_queryw <- stats::fft(data_fft_queryw_fft, inverse = TRUE) / length(data_fft_queryw_fft)

  last_product <- data_queryw[window_size:data_size]

  distance_profile <- data_pre - 2 * (last_product - sumwy * data_mean) / data_sd + sumwy2

  return(list(distance_profile = distance_profile, last_product = last_product))
}
