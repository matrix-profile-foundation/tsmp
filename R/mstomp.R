#' Multivariate STOMP algorithm
#'
#' No description
#'
#' No details
#'
#' @param data a matrix, where each colums is a time series (dimention). Can accept Lists and data.frames too.
#' @param query.size size of the sliding window
#' @param must.dim which dimentions to forcibly include (default is NULL)
#' @param exc.dim which dimentions to exclude (default is NULL)
#'
#' @return The matrix profile and profile index
#' @export
#'
#' @seealso [stamp()]
#' @references Yeh CM, Kavantzas N, Keogh E. Matrix Profile VI : Meaningful Multidimensional Motif Discovery.
#'
#' @examples
#' \dontrun{
#' mp <- mstomp(data, 30)
#' mp <- mstomp(data, 30, must.dim = c(1, 2))
#' mp <- mstomp(data, 30, exc.dim = c(2,3))
#' }
mstomp <- function(data, query.size, must.dim = NULL, exc.dim = NULL) {
  ## get various length
  exc_zone <- round(sub_len / 2)

  ## transform data list into matrix
  if (is.list(data)) {
    data_len <- length(data[[1]])
    n_dim <- length(data)

    for (i in 1:n_dim) {
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
    data_len <- nrow(data)
    n_dim <- ncol(data)
  } else if (is.vector(data)) {
    data_len <- length(data)
    n_dim <- 1
    # transform data into 1-col matrix
    data <- as.matrix(data) # just to be uniform
  } else {
    stop("Unknown type of data. Must be: matrix, data.frame, vector or list")
  }

  pro_len <- data_len - sub_len + 1

  ## check input
  if (sub_len > data_len / 2) {
    stop("Error: Time series is too short relative to desired subsequence length")
  }
  if (sub_len < 4) {
    stop("Error: Subsequence length must be at least 4")
  }
  if (any(must_dim > n_dim)) {
    stop("Error: The must have dimension must be less then the total dimension")
  }
  if (any(exc_dim > n_dim)) {
    stop("Error: The exclusion dimension must be less then the total dimension")
  }
  if (length(intersect(must_dim, exc_dim)) > 0) {
    stop("Error: The same dimension is presented in both the exclusion dimension and must have dimension")
  }

  ## check skip position
  n_exc <- length(exc_dim)
  n_must <- length(must_dim)
  mask_exc <- rep(FALSE, n_dim)
  mask_exc[exc_dim] <- TRUE
  skip_loc <- rep(FALSE, pro_len)

  for (i in 1:pro_len) {
    if (any(is.na(data[i:(i + sub_len - 1), !mask_exc])) || any(is.infinite(data[i:(i + sub_len - 1), !mask_exc]))) {
      skip_loc[i] <- TRUE
    }
  }

  data[is.na(data)] <- 0
  data[is.infinite(data)] <- 0

  pb <- txtProgressBar(min = 0, max = pro_len, style = 3)
  on.exit(close(pb))
  on.exit(beep(), TRUE)

  ## initialization
  data_freq <- matrix(0, (sub_len + data_len), n_dim)
  data_mu <- matrix(0, pro_len, n_dim)
  data_sig <- matrix(0, pro_len, n_dim)
  first_prod <- matrix(0, pro_len, n_dim)

  for (i in 1:n_dim) {
    nnPre <- mass_pre_mstamp(data[, i], data_len, sub_len)
    data_freq[, i] <- nnPre$data_freq
    data_mu[, i] <- nnPre$data_mu
    data_sig[, i] <- nnPre$data_sig

    mstomp <- mass_mstamp(data_freq[, i], data[1:sub_len, i], data_len, sub_len, data_mu[, i], data_sig[, i], data_mu[1, i], data_sig[1, i])
    first_prod[, i] <- mstomp$last_prod
  }

  ## compute the matrix profile
  pro_mul <- matrix(0, pro_len, n_dim)
  pro_idx <- matrix(0, pro_len, n_dim)
  dist_pro <- matrix(0, pro_len, n_dim)
  last_prod <- matrix(0, pro_len, n_dim)
  drop_val <- matrix(0, 1, n_dim)
  for (i in 1:pro_len) {
    # compute the distance profile
    setTxtProgressBar(pb, i)

    query <- as.matrix(data[i:(i + sub_len - 1), ])

    if (i == 1) {
      for (j in 1:n_dim) {
        mstomp <- mass_mstamp(data_freq[, j], query[, j], data_len, sub_len, data_mu[, j], data_sig[, j], data_mu[i, j], data_sig[i, j])
        dist_pro[, j] <- mstomp$dist_pro
        last_prod[, j] <- mstomp$last_prod
      }
    } else {

      rep_drop_val <- kronecker(matrix(1, pro_len - 1, 1), t(drop_val))
      rep_query <- kronecker(matrix(1, pro_len - 1, 1), t(query[sub_len, ]))

      last_prod[2:(data_len - sub_len + 1), ] <- last_prod[1:(data_len - sub_len), ] -
        data[1:(data_len - sub_len), ] * rep_drop_val +
        data[(sub_len + 1):data_len, ] * rep_query


      last_prod[1, ] <- first_prod[i, ]

      dist_pro <- 2 * (sub_len - (last_prod - sub_len * data_mu * kronecker(matrix(1, pro_len, 1), t(data_mu[i, ]))) /
                         (data_sig * kronecker(matrix(1, pro_len, 1), t(data_sig[i, ]))))
    }

    dist_pro <- Re(dist_pro)
    drop_val <- query[1, ]

    # apply exclusion zone
    exc_st <- max(1, i - exc_zone)
    exc_ed <- min(pro_len, i + exc_zone)
    dist_pro[exc_st:exc_ed, ] <- Inf
    dist_pro[data_sig < eps] <- Inf
    if (skip_loc[i] || any(data_sig[i, !mask_exc] < eps)) {
      dist_pro[] <- Inf
    }
    dist_pro[skip_loc, ] <- Inf

    # apply dimension "must have" and "exclusion"
    dist_pro[, exc_dim] <- Inf

    if (n_must > 0) {
      mask_must <- rep(FALSE, n_must)
      mask_must[must_dim] <- TRUE
      dist_pro_must <- dist_pro[, mask_must];
      dist_pro[, mask_must] <- -Inf
    }

    if (n_dim > 1)
      dist_pro_sort <- t(apply(dist_pro, 1, sort)) # sort by row, put all -Inf to the first columns
    else
      dist_pro_sort <- dist_pro

    if (n_must > 0)
      dist_pro_sort[, 1:n_must] <- dist_pro_must

    # figure out and store the nearest neighbor
    dist_pro_cum <- rep(0, pro_len)
    dist_pro_merg <- rep(0, pro_len)
    for (j in (max(1, n_must):(n_dim - n_exc))) {
      dist_pro_cum <- dist_pro_cum + dist_pro_sort[, j]
      dist_pro_merg[] <- dist_pro_cum / j
      min_idx <- which.min(dist_pro_merg)
      min_val <- dist_pro_merg[min_idx]
      pro_mul[i, j] <- min_val
      pro_idx[i, j] <- min_idx
    }
  }

  ## remove bad k setting in the returned matrix
  pro_mul <- sqrt(pro_mul)
  if (n_must > 0)
    pro_mul[, 1:(n_must - 1)] <- NA
  if (n_exc > 0)
    pro_mul[, (n_dim - n_exc + 1):length(pro_mul)] <- NA
  if (n_must > 0)
    pro_idx[, 1:(n_must - 1)] <- NA
  if (n_exc > 0)
    pro_idx[, (n_dim - n_exc + 1):length(pro_idx)] <- NA

  return(list(pro_mul = pro_mul, pro_idx = pro_idx))
}
