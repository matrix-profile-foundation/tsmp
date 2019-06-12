#' MP Distance
#'
#' @param ref_data
#' @param query_data
#' @param window_size
#'
#' @return
#' @export
#'
#' @examples
mpdist <- function(ref_data, query_data, window_size) {
  # TODO: check bugs in mass_v3 for large window_sizes
  if (length(ref_data) == length(query_data)) {
    return(mpdist_simple(ref_data, query_data, window_size))
  } else {
    return(mpdist_vect(ref_data, query_data, window_size))
  }
}

#' MP Distance for data of equal sizes
#'
#' @param ref_data reference data
#' @param query_data query data
#' @param window_size window size
#'
#' @return Returns the distance
#' @examples
#' dist <- mpdist(mp_toy_data$data[, 1], mp_toy_data$data[, 2], mp_toy_data$sub_len)
mpdist_simple <- function(ref_data, query_data, window_size) {
  thr <- 0.05
  mp <- mpx_abba_stomp(ref_data, query_data, window_size = window_size)

  dist <- cal_mp_dist(c(mp$mpa, mp$mpb), thr, length(ref_data) + length(query_data))

  return(dist)
}

#' MP Distance for data of different sizes
#'
#' @param data long data
#' @param query Short data
#' @param window_size window size
#'
#' @return Returns the distances for all data
#' @examples
mpdist_vect <- function(data, query, window_size) {
  thr <- 0.05
  mat <- NULL
  nn <- NULL

  num_subseqs <- length(query) - window_size + 1
  dist_profile_size <- length(data) - window_size + 1
  mat <- matrix(nrow = num_subseqs, ncol = dist_profile_size)

  for (i in seq_len(num_subseqs)) {
    nn <- dist_profile(data, query, nn, window_size = window_size, index = i, method = "v2")
    mat[i, ] <- Re(sqrt(nn$distance_profile))
  }

  all_right_histogram <- do.call(pmin, as.data.frame(t(mat))) # col min

  mass_dist_slid_min <- matrix(nrow = nrow(mat), ncol = ncol(mat))

  for (i in seq_len(nrow(mat))) {
    mass_dist_slid_min[i, ] <- caTools::runmin(mat[i, ], nrow(mat))
  }

  mp_dist_length <- length(data) - length(query) + 1
  right_hist_length <- length(query) - window_size + 1
  mpdist_array <- NULL # mp_dist_length - 1
  left_hist <- NULL # rightHistLength


  mpdist_array <- vector(mode = "numeric", length = mp_dist_length)
  for (i in seq_len(mp_dist_length)) {
    right_hist <- all_right_histogram[i:(right_hist_length + i - 1)]
    left_hist <- mass_dist_slid_min[, (i + floor(nrow(mat) / 2))]
    recreated_mp <- c(left_hist, right_hist)
    mpdist_array[i] <- cal_mp_dist(recreated_mp, thr, 2 * length(query))
  }

  return(mpdist_array)
}

#' Calculates the distance
#'
#' @param mp a matrix profile
#' @param thr threshold
#' @param data_size size of the data
#'
#' @return Returns the distance
#'
#' @examples
cal_mp_dist <- function(mp, thr, data_size, partial = TRUE) {
  k <- ceiling(thr * data_size)

  if (k > length(mp)) {
    # last element of a sorted list is the max()
    return(max(mp))
  }

  # keep classic sorting for debug reasons
  if (partial) {
    mp_sorted <- sort.int(mp, partial = k)
  } else {
    mp_sorted <- sort.int(mp)
  }

  dist <- mp_sorted[k]

  return(dist)
}

#' ABBA Stomp
#'
#' @param data
#' @param query
#' @param window_size
#'
#' @return
#'
#' @examples
mpx_abba_stomp <- function(data, query, window_size) {
  # forward
  nn <- dist_profile(data, query, window_size = window_size, method = "v2")
  # reverse
  # This is needed to handle with the join similarity.
  rnn <- dist_profile(query, data, window_size = window_size, method = "v2")

  data_size <- length(data)
  query_size <- length(query)
  matrix_profile_size <- data_size - window_size + 1
  matrix_profile2_size <- query_size - window_size + 1
  matrix_profile <- rep(Inf, matrix_profile_size)
  matrix_profile2 <- rep(Inf, matrix_profile2_size)
  num_queries <- query_size - window_size + 1

  first_product <- rnn$last_product

  for (i in 1:num_queries) {
    # compute the distance profile
    query_window <- query[i:(i + window_size - 1)]

    if (i == 1) {
      distance_profile <- nn$distance_profile
      last_product <- nn$last_product
    } else {
      last_product[2:(data_size - window_size + 1)] <- last_product[1:(data_size - window_size)] -
        data[1:(data_size - window_size)] * drop_value +
        data[(window_size + 1):data_size] * query_window[window_size]

      last_product[1] <- first_product[i]
      distance_profile <- 2 * (window_size - (last_product - window_size * nn$par$data_mean * nn$par$query_mean[i]) /
        (nn$par$data_sd * nn$par$query_sd[i]))
    }

    distance_profile <- Re(sqrt(distance_profile))
    drop_value <- query_window[1]

    ind <- (distance_profile < matrix_profile)
    matrix_profile[ind] <- distance_profile[ind]
    matrix_profile2[i] <- min(distance_profile)
  }

  return(list(mpa = matrix_profile, mpb = matrix_profile2))
}
