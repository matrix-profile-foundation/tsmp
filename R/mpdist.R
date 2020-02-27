#' MPdist - Distance between Time Series using Matrix Profile
#'
#' MPdist is a recently introduced distance measure which considers two time series to be similar
#' if they share many similar subsequences, regardless of the order of matching subsequences.
#' It was demonstrated in that MPdist is robust to spikes, warping, linear trends, dropouts,
#' wandering baseline and missing values, issues that are common outside of benchmark datasets.
#'
#' @details MPdist returns the distance of two time series or a vector containing the distance
#' between all sliding windows. If argument `type` is set to `vector`, the vector is returned.
#'
#' @param ref_data a `matrix` or a `vector`. The reference data
#' @param query_data a `matrix` or a `vector`. The query data
#' @param window_size an int. Size of the sliding window.
#' @param type the type of result. (Default is `simple`). See details.
#' @param thr threshold for MPdist. (Default is `0.05`). Don't change this unless you know what you are doing.
#'
#' @return Returns the distance of two time series or a vector containing the distance
#' between all sliding windows.
#' @export
#'
#' @references * Gharghabi S, Imani S, Bagnall A, Darvishzadeh A, Keogh E. Matrix Profile XII:
#' MPdist: A Novel Time Series Distance Measure to Allow Data Mining in More Challenging Scenarios.
#' In: 2018 IEEE International Conference on Data Mining (ICDM). 2018.
#'
#' @references Website: <https://sites.google.com/site/mpdistinfo/>
#'
#' @family distance measure
#'
#' @examples
#' ref_data <- mp_toy_data$data[, 1]
#' qe_data <- mp_toy_data$data[, 2]
#' qd_data <- mp_toy_data$data[150:200, 1]
#' w <- mp_toy_data$sub_len
#'
#' # distance between data of same size
#' deq <- mpdist(ref_data, qe_data, w)
#'
#' # distance between data of different sizes
#' ddiff <- mpdist(ref_data, qd_data, w)
#'
#' # distance vector between data of different sizes
#' ddvect <- mpdist(ref_data, qd_data, w, type = "vector")
mpdist <- function(ref_data, query_data, window_size, type = c("simple", "vector"), thr = 0.05) {
  type <- match.arg(type)

  # transform data into matrix
  if (is.vector(ref_data)) {
    ref_data <- as.matrix(ref_data)
  }
  else if (is.matrix(ref_data)) {
    if (ncol(ref_data) > nrow(ref_data)) {
      ref_data <- t(ref_data)
    }
  } else {
    stop("Unknown type of data. Must be: a column matrix or a vector.", call. = FALSE)
  }

  if (is.vector(query_data)) {
    query_data <- as.matrix(query_data)
  } else if (is.matrix(query_data)) {
    if (ncol(query_data) > nrow(query_data)) {
      query_data <- t(query_data)
    }
  } else {
    stop("Unknown type of query. Must be: a column matrix or a vector.", call. = FALSE)
  }

  if (window_size < 4) {
    stop("`window_size` must be at least 4.", call. = FALSE)
  }

  if (nrow(ref_data) < nrow(query_data)) {
    temp <- ref_data
    ref_data <- query_data
    query_data <- temp
  }

  if (type == "simple" || nrow(ref_data) == nrow(query_data)) {
    if (nrow(query_data) == window_size) {
      warning("Distance profile is being used since window_size equals to query size")
      dist <- dist_profile(ref_data, query_data)$distance_profile
      dist <- sqrt(cal_mp_dist(dist, thr, nrow(ref_data)))
    } else {
      dist <- mpdist_simple(ref_data, query_data, window_size, thr)
    }
  } else {
    obj <- list()
    if (nrow(query_data) == window_size) {
      warning("Distance profile is being used since window_size equals to query size")
      obj$mpdist <- sqrt(dist_profile(ref_data, query_data)$distance_profile)
    } else {
      obj$mpdist <- mpdist_vect(ref_data, query_data, window_size, thr)
    }
    obj$w <- window_size
    obj$data <- list(ref_data, query_data)
    class(obj) <- "MPdistProfile"

    attr(obj, "origin") <- list(
      data_size = nrow(ref_data),
      query_size = nrow(query_data),
      window_size = window_size,
      mp_size = nrow(obj$mpdist),
      algorithm = "MPdist",
      class = class(obj),
      version = 1.1
    )

    return(obj)
  }

  return(dist)
}

#' MP Distance
#'
#' @param data reference data
#' @param query query data
#' @param window_size window size
#' @param thr threshold for mpdist
#'
#' @return Returns the distance vector
#'
#' @keywords internal
#' @noRd
mpdist_simple <- function(ref_data, query_data, window_size, thr = 0.05) {
  mp <- mpx(data = ref_data, query = query_data, window_size = window_size, idx = FALSE)

  dist <- cal_mp_dist(c(mp$mpa, mp$mpb), thr, nrow(ref_data) + nrow(query_data))

  return(dist)
}

#' MP Distance vector
#'
#' @param data reference data
#' @param query query data
#' @param window_size window size
#' @param thr threshold for mpdist
#'
#' @return Returns the distance vector
#'
#' @keywords internal
#' @noRd
mpdist_vect <- function(data, query, window_size, thr = 0.05) {
  nn <- NULL
  query_size <- nrow(query)
  data_size <- nrow(data)
  num_subseqs <- query_size - window_size + 1
  dist_profile_size <- data_size - window_size + 1
  mat <- matrix(nrow = num_subseqs, ncol = dist_profile_size)

  for (i in seq_len(num_subseqs)) {
    nn <- dist_profile(data, query, nn, window_size = window_size, index = i)
    mat[i, ] <- nn$distance_profile
  }

  all_right_histogram <- do.call(pmin, as.data.frame(t(mat))) # col min

  mass_dist_slid_min <- matrix(nrow = num_subseqs, ncol = dist_profile_size - num_subseqs + 1)
  # apply is not faster
  for (i in seq_len(num_subseqs)) {
    mass_dist_slid_min[i, ] <- movmin(mat[i, ], num_subseqs)
  }

  mp_dist_length <- data_size - query_size + 1
  right_hist_length <- query_size - window_size + 1 # num_subseqs
  mpdist_array <- vector(mode = "numeric", length = mp_dist_length)

  for (i in seq_len(mp_dist_length)) {
    right_hist <- all_right_histogram[i:(right_hist_length + i - 1)]
    left_hist <- mass_dist_slid_min[, i]
    recreated_mp <- c(left_hist, right_hist)
    mpdist_array[i] <- cal_mp_dist(recreated_mp, thr, 2 * query_size)
  }

  mpdist_array[mpdist_array < vars()$eps] <- 0
  mpdist_array <- sqrt(mpdist_array)

  obj <- as.matrix(mpdist_array)

  return(obj)
}

#' Calculates the distance
#'
#' @param mp a matrix profile
#' @param thr threshold for mpdist
#' @param data_size size of the data
#'
#' @return Returns the distance
#'
#' @keywords internal
#' @noRd
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
