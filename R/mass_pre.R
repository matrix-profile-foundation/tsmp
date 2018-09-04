#' Precomputes several values used on MASS
#'
#' @param data a `vector` or a `matrix` of `numeric`. Reference Time Series.
#' @param data.size an `int`. Reference Time Series size.
#' @param query a `vector` or a `matrix` of `numeric`. Query Time Series (default is `NULL`).
#' @param query.size an `int`. Query Time Series size (default is `NULL`).
#' @param window.size an `int`. Sliding window size.
#'
#' @return Returns `data.fft`, `data.mean`, `data.sd`, `query.mean` and `query.sd`.
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
