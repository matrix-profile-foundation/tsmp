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

mass.pre <- function(data, data.size, query = NULL, query.size = NULL, window.size) {

  if (is.matrix(data))
    data <- as.vector(data)

  data.mean <- fast.movavg(data, window.size) # precompute moving average
  data.sd <- fast.movsd(data, window.size) # precompute moving SD
  data[(data.size + 1):(window.size + data.size)] <- 0
  data.fft <- stats::fft(data) # precompute fft of data


  if (!is.null(query)) {
    if (is.matrix(query))
      query <- as.vector(query)

    query.mean <- fast.movavg(query, window.size) # precompute moving average
    query.sd <- fast.movsd(query, window.size) # precompute moving SD
  } else {
    query.mean <- data.mean
    query.sd <- data.sd
  }


  return(list(data.fft = data.fft, data.mean = data.mean, data.sd = data.sd, query.mean = query.mean, query.sd = query.sd))
}
