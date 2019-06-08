#' Title
#'
#' @param data
#' @param query
#' @param window_size
#'
#' @return
#' @export
#'
#' @examples
#'

mpdist <- function(data, query, window_size) {
  thr <- 0.05
  mpa <- stomp(data, query, window_size = window_size, verbose = 0) # TODO: mpx_ABBA(a, b, w)
  mpb <- stomp(query, data, window_size = window_size, verbose = 0)

  dist <- cal_mp_dist(c(mpa$mp, mpb$mp), thr, 2 * length(query))

  return(dist)
}

mpdist_vect <- function(data, query, window_size) {
  vec <- NULL
  for (i in seq_len(length(data) - length(query) + 1)) {
    vec <- c(vec, mpdist(query, data[i:(i + length(query) - 1)], window_size))
  }

  return(vec)
}

fast_mpdist_vect <- function(data, query, window_size) {
  thr <- 0.05
  mat <- NULL
  nn <- NULL

  for (i in seq_len(length(query) - window_size + 1)) {
    nn <- dist_profile(data, query, nn, window_size = window_size, index = i, method = "v3")
    mat <- rbind(mat, Re(sqrt(nn$distance_profile)))
  }

  all_right_histogram <- do.call(pmin, as.data.frame(t(mat))) # col min


  mass_dist_slid_min <- NULL

  for (i in seq_len(nrow(mat))) {
    mass_dist_slid_min <- rbind(
      mass_dist_slid_min,
      caTools::runmin(mat[i, ], window_size)
    )
  }

  mp_dist_length <- length(data) - length(query) + 1
  right_hist_length <- length(query) - window_size + 1
  mpdist_array <- NULL # mp_dist_length - 1
  left_hist <- NULL # rightHistLength


  for (i in seq_len(mp_dist_length)) {
    right_hist <- all_right_histogram[i:(right_hist_length + i - 1)]
    left_hist <- mass_dist_slid_min[, (i + floor(nrow(mat) / 2))]
    recreated_mp <- c(left_hist, right_hist)
    mpdist_array <- c(mpdist_array, cal_mp_dist(recreated_mp, thr, 2 * length(query)))
  }

  return(mpdist_array)
}

cal_mp_dist <- function(matrix_profile, thr, data_size) {
  mp_sorted <- sort(matrix_profile)
  k <- ceiling(thr * data_size)

  if (k <= length(mp_sorted)) {
    dist <- mp_sorted[k]
  } else {
    dist <- tail(mp_sorted, 1)
  }

  return(dist)
}
