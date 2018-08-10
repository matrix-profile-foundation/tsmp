#' Title
#'
#' @param data
#' @param sub_len
#' @param pro_mul
#' @param pro_idx
#' @param n_dim
#'
#' @return
#' @export
#'
#' @examples
guide_search <- function(data, sub_len, pro_mul, pro_idx, n_dim) {

  ## transform data list into matrix
  if (is.list(data)) {
    data_len <- length(data[[1]])
    data_dim <- length(data)

    for (i in 1:data_dim) {
      len <- length(data[[i]])
      # Fix TS size with NaN
      if (len < data_len) {
        data[[i]] <- c(data[[i]], rep(NA, data_len - len))
      }
    }
    # transform data into matrix (each column is a TS)
    data <- sapply(data, cbind)
  } else if (is.matrix(data) || is.data.frame(data)) {
    if (is.data.frame(data)) {
      data <- as.matrix(data)
    } # just to be uniform
    if (ncol(data) > nrow(data))
      data <- t(data)
    data_len <- nrow(data)
    n_dim <- ncol(data)
  } else if (is.vector(data)) {
    data_len <- length(data)
    # transform data into 1-col matrix
    data <- as.matrix(data) # just to be uniform
  } else {
    stop("Unknown type of data. Must be: matrix, data.frame, vector or list")
  }

  pro_mul <- pro_mul[, n_dim]
  pro_idx <- pro_idx[, n_dim]
  motif_idx <- which.min(pro_mul)
  motif_idx <- sort(c(motif_idx, pro_idx[motif_idx]))

  motif_1 <- data[motif_idx[1]:(motif_idx[1] + sub_len - 1), ]
  motif_2 <- data[motif_idx[2]:(motif_idx[2] + sub_len - 1), ]

  motif_dim <- sort(sum(abs(motif_1 - motif_2), 1), index.return = TRUE)$ix
  motif_dim <- sort(motif_dim[1:n_dim])
  motif_dim <- list(motif_dim, motif_dim)

  return(list(motif_idx = motif_idx, motif_dim = motif_dim))
}
