#' Precomputes several values used on MASS
#'
#' @param data a `vector` or a `matrix` of `numeric`. Reference Time Series.
#' @param query a `vector` or a `matrix` of `numeric`. Query Time Series (default is `NULL`).
#' @param window_size an `int`. Sliding window size.
#' @param weight a `vector` of `numeric` with the same length of the `window_size`.
#'
#' @return Returns `window_size`, `data_fft`, `data_size`, `data_mean`, `data_sd`, `data_pre`, `weight`
#'
#' @keywords internal
#'
#' @seealso [mass_weighted()] for using precomputed values.
#'
#' @references * Abdullah Mueen, Yan Zhu, Michael Yeh, Kaveh Kamgar, Krishnamurthy Viswanathan,
#'   Chetan Kumar Gupta and Eamonn Keogh (2015), The Fastest Similarity Search Algorithm for Time
#'   Series Subsequences under Euclidean Distance.
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
mass_pre_w <- function(data, query = NULL, window_size, weight) {
  if (is.matrix(data)) {
    data <- as.vector(data)
  }

  data_size <- length(data)

  if (window_size > data_size) {
    stop("'window_size' must be smaller or equal to 'data' length.")
  }

  if (is.matrix(weight)) {
    weight <- as.vector(weight)
  }

  if (length(weight) != window_size) {
    stop("'weight' must have the same length as the 'window_size'.")
  }

  data_avgsd <- fast_avg_sd(data, window_size) # precompute moving average and SD
  data_mean <- data_avgsd$avg
  data_sd <- data_avgsd$sd
  data_fft <- stats::fft(data) # precompute fft of data

  sumw <- sum(weight)
  rweight <- rev(weight)
  rweight[(window_size + 1):data_size] <- 0
  w_fft <- stats::fft(rweight)
  data_fft_w_fft <- data_fft * w_fft
  data_w <- stats::fft(data_fft_w_fft, inverse = TRUE) / length(data_fft_w_fft)
  data2_fft <- stats::fft(data^2)
  data2_fft_w_fft <- data2_fft * w_fft
  data2_w <- stats::fft(data2_fft_w_fft, inverse = TRUE) / length(data2_fft_w_fft)
  sumxw2 <- data2_w[window_size:data_size]
  sumxw <- data_w[window_size:data_size]

  data_pre <- ((sumxw2 - 2 * sumxw * data_mean + sumw * data_mean^2) / data_sd^2)



  if (!is.null(query)) {
    if (is.matrix(query)) {
      query <- as.vector(query)
    }

    if (window_size > length(query)) {
      stop("'window_size' must be smaller or equal to 'query' length.")
    }
  }

  return(list(
    window_size = window_size, data_fft = data_fft, data_size = data_size,
    data_mean = data_mean, data_sd = data_sd, data_pre = data_pre, weight = weight
  ))
}
