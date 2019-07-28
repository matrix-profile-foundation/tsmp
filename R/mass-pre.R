#' Precomputes several values used on MASS
#'
#' @param data a `vector` or a `matrix` of `numeric`. Reference Time Series.
#' @param query a `vector` or a `matrix` of `numeric`. Query Time Series (default is `NULL`).
#' @param window_size an `int`. Sliding window size.
#'
#' @return Returns `window_size`, `data_fft`, `data_size`, `data_mean`, `data_sd`, `query_mean` and `query_sd`.
#' @export
#' @keywords internal
#'
#' @seealso [mass_v2()], [mass_v3()] for using precomputed values.
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
#'
#' pre <- tsmp:::mass_pre(ref_data, query_data, w)
#'
#' dp <- list()
#' for (i in 1:(pre$data_size - w + 1)) {
#'   dp[[i]] <- tsmp:::mass_v2(
#'     query_data[i:(i - 1 + w)], pre$window_size, pre$data_fft, pre$data_size,
#'     pre$data_mean, pre$data_sd, pre$query_mean[i], pre$query_sd[i]
#'   )
#' }
mass_pre <- function(data, query = NULL, window_size) {
  if (is.matrix(data)) {
    data <- as.vector(data)
  }

  data_size <- length(data)

  if (window_size > data_size) {
    stop("'window_size' must be smaller or equal to 'data' length.")
  }

  data_avgsd <- fast_avg_sd(data, window_size) # precompute moving average and SD
  data_mean <- data_avgsd$avg
  data_sd <- data_avgsd$sd
  pad_size <- 2^ceiling(log2(data_size))
  data[(data_size + 1):pad_size] <- 0
  data_fft <- stats::fft(data) # precompute fft of data


  if (!is.null(query)) {
    if (is.matrix(query)) {
      query <- as.vector(query)
    }

    query_size <- length(query)

    if (window_size > query_size) {
      stop("'window_size' must be smaller or equal to 'query' length.")
    }

    query_avgsd <- fast_avg_sd(query, window_size) # precompute moving average and SD
    query_mean <- query_avgsd$avg
    query_sd <- query_avgsd$sd
  } else {
    query_mean <- data_mean
    query_sd <- data_sd
  }

  return(list(
    window_size = window_size, data_fft = data_fft, data_size = data_size,
    data_mean = data_mean, data_sd = data_sd, query_mean = query_mean,
    query_sd = query_sd
  ))
}
