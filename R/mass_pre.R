#' Precomputes several values used on MASS
#'
#' @param data reference Time Series
#' @param data.size reference Time Series size
#' @param query query Time Series (default is NULL)
#' @param query.size query Time Series (default is NULL)
#' @param window.size sliding window size
#'
#' @return Returns data.fft, data.mean, data.sd, query.mean and query.sd
#' @export
#'
#' @references Abdullah Mueen, Yan Zhu, Michael Yeh, Kaveh Kamgar, Krishnamurthy Viswanathan, Chetan Kumar Gupta and Eamonn Keogh (2015), The Fastest Similarity Search Algorithm for Time Series Subsequences under Euclidean Distance
#' @references <https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html>

mass.pre <- function(data, data.size, query = NULL, query.size = NULL, window.size) {

  if (is.matrix(data))
    data <- as.vector(data)

  data.mean <- fast.movavg(data, window.size) # precompute moving average
  data.sd <- fast.movsd(data, window.size) # precompute moving SD
  data[(data.size + 1):(window.size + data.size)] <- 0
  data.fft <- fft(data) # precompute fft of data


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
