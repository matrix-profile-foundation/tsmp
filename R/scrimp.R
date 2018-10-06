#' Anytime univariate SCRIMP algorithm (experimental)
#'
#' Computes the best so far Matrix Profile and Profile Index for Univariate Time Series.
#' DISCLAIMER: This algorithm still in development by its authors.
#' Join similarity, RMP and LMP not implemented yet.
#'
#' @details
#' The Matrix Profile, has the potential to revolutionize time series data mining because of its
#' generality, versatility, simplicity and scalability. In particular it has implications for time
#' series motif discovery, time series joins, shapelet discovery (classification), density
#' estimation, semantic segmentation, visualization, rule discovery, clustering etc. The anytime
#' SCRIMP computes the Matrix Profile and Profile Index in such manner that it can be stopped before
#' its complete calculation and return the best so far results allowing ultra-fast approximate
#' solutions. `verbose` changes how much information is printed by this function; `0` means nothing,
#' `1` means text, `2` adds the progress bar, `3` adds the finish sound. `exclusion_zone` is used to
#' avoid  trivial matches.
#'
#' @param ... a `matrix` or a `vector`.
#' @param window_size an `int`. Size of the sliding window.
#' @param exclusion_zone a `numeric`. Size of the exclusion zone, based on window size (default is
#'   `1/2`). See details.
#' @param verbose an `int`. See details. (Default is `2`).
#' @param s_size a `numeric`. for anytime algorithm, represents the size (in observations) the
#'   random calculation will occur (default is `Inf`).
#'
#' @return Returns a `MatrixProfile` object, a `list` with the matrix profile `mp`, profile index `pi`
#'   left and right matrix profile `lmp`, `rmp` and profile index `lpi`, `rpi`, window size `w` and
#'   exclusion zone `ez`.
#' @export
#'
#' @family matrix profile computations
#'
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#' mp <- scrimp(mp_toy_data$data[1:200, 1], window_size = 30, verbose = 0)
#'
#' \dontrun{
#' ref_data <- mp_toy_data$data[, 1]
#' query_data <- mp_toy_data$data[, 2]
#' # self similarity
#' mp <- scrimp(ref_data, window_size = 30, s_size = round(nrow(ref_data) * 0.1))
#' # join similarity
#' mp <- scrimp(ref_data, query_data, window_size = 30, s_size = round(nrow(query_data) * 0.1))
#' }

scrimp <- function(..., window_size, exclusion_zone = 1 / 2, verbose = 2, s_size = Inf) {
  args <- list(...)
  data <- args[[1]]
  if (length(args) > 1) {
    message("Join similarity not implemented yet.")
    # invisible(return(NULL))
    query <- args[[2]]
    exclusion_zone <- 0 # don't use exclusion zone for joins
  } else {
    query <- data
  }

  # transform data into matrix
  if (is.vector(data)) {
    data <- as.matrix(data)
  }
  else if (is.matrix(data)) {
    if (ncol(data) > nrow(data)) {
      data <- t(data)
    }
  } else {
    stop("Error: Unknown type of data. Must be: a column matrix or a vector.")
  }

  if (is.vector(query)) {
    query <- as.matrix(query)
  } else if (is.matrix(query)) {
    if (ncol(query) > nrow(query)) {
      query <- t(query)
    }
  } else {
    stop("Error: Unknown type of query. Must be: a column matrix or a vector.")
  }

  ez <- exclusion_zone # store original
  exclusion_zone <- round(window_size * exclusion_zone + vars()$eps)
  data_size <- nrow(data)
  query_size <- nrow(query)

  if (query_size > data_size) {
    stop("Error: Query must be smaller or the same size as reference data.")
  }
  if (window_size > query_size / 2) {
    stop("Error: Time series is too short relative to desired window size.")
  }
  if (window_size < 4) {
    stop("Error: `window_size` must be at least 4.")
  }

  matrix_profile_size <- data_size - window_size + 1
  num_queries <- query_size - window_size + 1

  # check skip position
  skip_location <- rep(FALSE, matrix_profile_size)

  for (i in 1:matrix_profile_size) {
    if (any(is.na(data[i:(i + window_size - 1)])) || any(is.infinite(data[i:(i + window_size - 1)]))) {
      skip_location[i] <- TRUE
    }
  }

  data[is.na(data)] <- 0
  data[is.infinite(data)] <- 0

  query[is.na(query)] <- 0
  query[is.infinite(query)] <- 0

  matrix_profile <- matrix(Inf, matrix_profile_size, 1)
  profile_index <- matrix(-1, matrix_profile_size, 1)

  if (length(args) > 1) {
    # no RMP and LMP for joins
    left_matrix_profile <- right_matrix_profile <- NULL
    left_profile_index <- right_profile_index <- NULL
  } else {
    left_matrix_profile <- right_matrix_profile <- matrix_profile
    left_profile_index <- right_profile_index <- profile_index
  }

  j <- 1
  if (exclusion_zone > 0) {
    exclusion_zone <- exclusion_zone + 1
  }
  order <- (exclusion_zone + 1):num_queries
  ssize <- min(s_size, length(order))
  order <- sample(order, size = ssize)

  tictac <- Sys.time()

  if (verbose > 1) {
    pb <- utils::txtProgressBar(min = 1, max = ssize, style = 3, width = 80)
    on.exit(close(pb))
  }
  if (verbose > 2) {
    on.exit(beep(sounds[[1]]), TRUE)
  }
  # anytime must return the result always
  on.exit(return({
    if (length(args) == 1) {
      right_matrix_profile <- sqrt(abs(right_matrix_profile))
      left_matrix_profile <- sqrt(abs(left_matrix_profile))
    }
    obj <- list(
      mp = sqrt(abs(matrix_profile)), pi = profile_index,
      rmp = right_matrix_profile, rpi = right_profile_index,
      lmp = left_matrix_profile, lpi = left_profile_index,
      w = window_size,
      ez = ez
    )
    class(obj) <- "MatrixProfile"
    obj
  }), TRUE)

  pre_data <- fast_avg_sd(data, window_size)
  pre_query <- fast_avg_sd(query, window_size)

  curlastz <- rep(0, num_queries)
  curdistance <- rep(0, num_queries)
  dist1 <- rep(Inf, num_queries)
  dist2 <- rep(Inf, num_queries)
  index_order <- seq_len(num_queries)

  for (i in order) {
    j <- j + 1

    curlastz[i] <- sum(data[1:window_size] * query[i:(i + window_size - 1)])

    curlastz[(i + 1):num_queries] <-
      curlastz[i] +
      cumsum(
        query[(i + window_size):data_size] * data[(window_size + 1):(query_size - i + 1)] # a_term
        - data[1:(num_queries - i)] * query[i:(num_queries - 1)] # m_term
      )

    curdistance[i:num_queries] <-
      2 * (window_size -
        (curlastz[i:num_queries] # x_term
        - window_size * pre_query$avg[i:num_queries] * pre_data$avg[1:(num_queries - i + 1)]) /
          (pre_query$sd[i:num_queries] * pre_data$sd[1:(num_queries - i + 1)])
      )

    # Skip positions
    skipped_curdistance <- curdistance
    # skipped_curdistance[pre_data$sd[i:num_queries] < vars()$eps] <- Inf
    # if (skip_location[i] || any(pre_query$sd[i] < vars()$eps)) {
    #   skipped_curdistance[] <- Inf
    # }
    # skipped_curdistance[skip_location[i:num_queries]] <- Inf

    # update matrix profile
    dist1[1:(i - 1)] <- Inf
    dist1[i:num_queries] <- skipped_curdistance[i:num_queries]
    dist2[1:(num_queries - i + 1)] <- skipped_curdistance[i:num_queries]
    dist2[(num_queries - i + 2):num_queries] <- Inf

    loc1 <- (dist1 < matrix_profile)
    matrix_profile[loc1] <- dist1[loc1]
    profile_index[loc1] <- index_order[loc1] - i + 1
    loc2 <- (dist2 < matrix_profile)
    matrix_profile[loc2] <- dist2[loc2]
    profile_index[loc2] <- index_order[loc2] + i - 1

    if (length(args) == 1) {
      # left matrix_profile
      loc1 <- (dist1 < left_matrix_profile)
      left_matrix_profile[loc1] <- dist1[loc1]
      left_profile_index[loc1] <- index_order[loc1] - i + 1

      # right matrix_profile
      loc2 <- (dist2 < right_matrix_profile)
      right_matrix_profile[loc2] <- dist2[loc2]
      right_profile_index[loc2] <- index_order[loc2] + i - 1
    }

    if (verbose > 1) {
      utils::setTxtProgressBar(pb, j)
    }
  }

  tictac <- Sys.time() - tictac

  if (verbose > 0) {
    message(sprintf("\nFinished in %.2f %s", tictac, units(tictac)))
  }

  # return() is at on.exit() function
}

#' SCRIMP NN algorithm
#'
#' @param data a `matrix` or a `vector`.
#' @param idx an `int`. The index of window to be computed.
#' @param data_size an `int`.
#' @param window_size an `int`.
#' @param mp_size an `int`.
#' @param data_mean result of `fast_movavg`.
#' @param data_sd result of `fast_movsd`.
#'
#' @return Returns the distance profile
#' @keywords internal
#' @noRd

diagonal_dist <- function(data, query, idx, data_size, query_size, window_size, mp_size, data_mean, data_sd, query_mean, query_sd) {
  data <- as.matrix(data)
  query <- as.matrix(query)
  x_term <- matrix(1, mp_size - idx + 1, 1) * (t(data[idx:(idx + window_size - 1), 1, drop = FALSE]) %*% query[1:window_size, 1, drop = FALSE])[1]
  m_term <- data[idx:(mp_size - 1)] * query[1:(mp_size - idx)]
  a_term <- data[(idx + window_size):data_size] * query[(window_size + 1):(query_size - idx + 1)]
  if (mp_size != idx) {
    x_term[2:length(x_term)] <- x_term[2:length(x_term)] - cumsum(m_term) + cumsum(a_term)
  }

  distance_profile <- (x_term - window_size * data_mean[idx:mp_size] * query_mean[1:(mp_size - idx + 1)]) /
    (window_size * data_sd[idx:length(data_sd)] * query_sd[1:(mp_size - idx + 1)])
  distance_profile <- 2 * window_size * (1 - distance_profile)

  return(distance_profile)
}
