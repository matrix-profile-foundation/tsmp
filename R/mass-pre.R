#' Precomputes several values used on MASS
#'
#' @param data a `vector` or a `matrix` of `numeric`. Reference Time Series.
#' @param data_size an `int`. Reference Time Series size.
#' @param query a `vector` or a `matrix` of `numeric`. Query Time Series (default is `NULL`).
#' @param query_size an `int`. Query Time Series size (default is `NULL`).
#' @param window_size an `int`. Sliding window size.
#'
#' @return Returns `data_fft`, `data_mean`, `data_sd`, `query_mean` and `query_sd`.
#' @export
#'
#' @seealso [mass()] for using precomputed values.
#'
#' @references * Abdullah Mueen, Yan Zhu, Michael Yeh, Kaveh Kamgar, Krishnamurthy Viswanathan,
#'   Chetan Kumar Gupta and Eamonn Keogh (2015), The Fastest Similarity Search Algorithm for Time
#'   Series Subsequences under Euclidean Distance.
#' @references Website: <https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html>
#'
#' @examples
#' w <- mp_toy_data$sub_len
#' ref_data <- mp_toy_data$data[,1]
#' query_data <- mp_toy_data$data[,1]
#' d_size <- length(ref_data)
#' q_size <- length(query_data)
#'
#' pre <- mass_pre(ref_data, d_size, query_data, q_size, w)
#'
#' dp <- list()
#' for(i in 1:(d_size - w + 1)) {
#'   dp[[i]] <- mass(pre$data_fft, query_data[i:(i-1+w)], d_size, w, pre$data_mean, pre$data_sd,
#'           pre$query_mean[i], pre$query_sd[i])
#' }

mass_pre <- function(data, data_size, query = NULL, query_size = NULL, window_size) {
  if (is.matrix(data)) {
    data <- as.vector(data)
  }

  data_mean <- fast_movavg(data, window_size) # precompute moving average
  data_sd <- fast_movsd(data, window_size) # precompute moving SD
  data[(data_size + 1):(window_size + data_size)] <- 0
  data_fft <- stats::fft(data) # precompute fft of data


  if (!is.null(query)) {
    if (is.matrix(query)) {
      query <- as.vector(query)
    }

    query_mean <- fast_movavg(query, window_size) # precompute moving average
    query_sd <- fast_movsd(query, window_size) # precompute moving SD
  } else {
    query_mean <- data_mean
    query_sd <- data_sd
  }


  return(list(data_fft = data_fft, data_mean = data_mean, data_sd = data_sd, query_mean = query_mean, query_sd = query_sd))
}
