#' Calculates the distance profile using MASS_V3 algorithm
#'
#' Mueen's Algorithm for Similarity Search is The Fastest Similarity Search Algorithm for Time
#' Series Subsequences under Euclidean Distance and Correlation Coefficient.
#'
#' This is a piecewise version of MASS that performs better when the size of the pieces are well
#' aligned with the hardware.
#'
#' @param data a `matrix` or a `vector`.
#' @param query_window a `vector` of `numeric`. Query window.
#' @param data_size an `int`. The length of the reference data.
#' @param window_size an `int`. Sliding window size.
#' @param data_mean precomputed data moving average.
#' @param data_sd precomputed data moving standard deviation.
#' @param query_mean precomputed query average.
#' @param query_sd precomputed query standard deviation.
#' @param k an `int` or `NULL`. Default is `NULL`. Defines the size of batch. Prefer to use a power of 2.
#' @param \dots just a placeholder to catch unused parameters.
#'
#' @return Returns the `distance_profile` for the given query and the `last_product` for STOMP
#'   algorithm.
#'
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
#' pre <- tsmp:::mass_pre(ref_data, query_data, w)
#'
#' dp <- list()
#' for (i in 1:(d_size - w + 1)) {
#'   dp[[i]] <- tsmp:::mass_v3(
#'     query_data[i:(i - 1 + w)], ref_data,
#'     pre$window_size, pre$data_size, pre$data_mean, pre$data_sd,
#'     pre$query_mean[i], pre$query_sd[i]
#'   )
#' }
mass_v3 <- function(query_window, data, window_size, data_size, data_mean, data_sd, query_mean, query_sd, k = NULL, ...) {
  distance_profile <- vector(mode = "numeric", data_size - window_size + 1)
  last_product <- vector(mode = "numeric", data_size - window_size + 1)

  if (is.null(k)) {
    k <- 1024
  }

  if (k > data_size) {
    k <- 2^ceiling(log2(sqrt(data_size)))
  }

  if (k <= window_size) {
    k <- 2^(ceiling(log2(window_size)) + 1)
    if (k > data_size) {
      k <- data_size
    }
  }

  # pre-process query for fft
  query_window <- rev(query_window)
  query_window[(window_size + 1):k] <- 0
  query_fft <- stats::fft(query_window)

  jump <- k - window_size + 1
  seq_end <- data_size - k + 1

  for (j in seq.int(1, seq_end, jump)) {
    idx_begin <- j
    idx_end <- j + k - window_size
    # The main trick of getting dot products in O(n log n) time
    data_fft <- stats::fft(data[j:(j + k - 1)])
    prod <- data_fft * query_fft
    z <- Re(stats::fft(prod, inverse = TRUE) / length(prod))
    d <- 2 * (window_size - (z[window_size:k] - window_size * data_mean[idx_begin:idx_end] * query_mean) / (data_sd[idx_begin:idx_end] * query_sd))
    distance_profile[idx_begin:idx_end] <- d
    last_product[idx_begin:idx_end] <- z[window_size:k]
  }

  j <- j + k - window_size + 1
  k <- data_size - j + 1


  if (k >= window_size) {
    idx_begin <- j
    idx_end <- data_size - window_size + 1
    # The main trick of getting dot products in O(n log n) time
    data_fft <- stats::fft(data[j:data_size])
    query_window <- query_window[1:k]
    query_fft <- stats::fft(query_window)
    prod <- data_fft * query_fft
    z <- Re(stats::fft(prod, inverse = TRUE) / length(prod))

    d <- 2 * (window_size - (z[window_size:k] - window_size * data_mean[idx_begin:idx_end] * query_mean) / (data_sd[idx_begin:idx_end] * query_sd))
    distance_profile[idx_begin:idx_end] <- d
    last_product[idx_begin:idx_end] <- z[window_size:k]
  }

  distance_profile[distance_profile < 0] <- 0

  return(list(distance_profile = distance_profile, last_product = last_product))
}
