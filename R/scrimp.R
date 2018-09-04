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
#' `1` means text, `2` means text and sound. `exclusion_zone` is used to avoid  trivial matches.
#'
#' @param ... a `matrix` or a `vector`.
#' @param window_size an `int`. Size of the sliding window.
#' @param exclusion_zone a `numeric`. Size of the exclusion zone, based on window size (default is
#'   `1/2`). See details.
#' @param s_size a `numeric`. for anytime algorithm, represents the size (in observations) the
#'   random calculation will occur (default is `Inf`).
#' @param verbose an `int`. See details. (Default is `2`).
#'
#' @return Returns the matrix profile `mp` and profile index `pi`.
#' @export
#'
#' @family matrix profile computations
#' @seealso [mstomp()], [mstomp_par()]
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#' mp <- scrimp(mp_toy_data$data[1:200,1], window_size = 30, verbose = 0)
#' \dontrun{
#' ref_data <- mp_toy_data$data[,1]
#' query_data <- mp_toy_data$data[,2]
#' # self similarity
#' mp <- scrimp(ref_data, window_size = 30, s_size = round(nrow(ref_data) * 0.1))
#' # join similarity
#' mp <- scrimp(ref_data, query_data, window_size = 30, s_size = round(nrow(query_data) * 0.1))
#' }

scrimp <- function(..., window_size, exclusion_zone = 1 / 2, s_size = Inf, verbose = 2) {
  args <- list(...)
  data <- args[[1]]
  if (length(args) > 1) {
    message("DISCLAIMER: This algorithm still in development by its authors.")
    message("Join similarity not implemented yet.")
    invisible(return(NULL))
    query <- args[[2]]
    exclusion_zone <- 0 # don't use exclusion zone for joins
  } else {
    query <- data
  }

  ## transform data into matrix
  if (is.vector(data)) {
    data <- as.matrix(data)
  }
  else if (is.matrix(data)) {
    if (ncol(data) > nrow(data)) {
      data <- t(data)
    }
  } else {
    stop("Error: Unknown type of data. Must be: a column matrix or a vector.", call. = FALSE)
  }

  if (is.vector(query)) {
    query <- as.matrix(query)
  } else if (is.matrix(query)) {
    if (ncol(query) > nrow(query)) {
      query <- t(query)
    }
  } else {
    stop("Error: Unknown type of query. Must be: a column matrix or a vector.", call. = FALSE)
  }

  exclusion_zone <- floor(window_size * exclusion_zone)
  data_size <- nrow(data)
  query_size <- nrow(query)
  matrix_profile_size <- data_size - window_size + 1
  num_queries <- query_size - window_size + 1

  ## check skip position
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

  if (query_size > data_size) {
    stop("Error: Query must be smaller or the same size as reference data.", call. = FALSE)
  }
  if (window_size > query_size / 2) {
    stop("Error: Time series is too short relative to desired window size.", call. = FALSE)
  }
  if (window_size < 4) {
    stop("Error: `window_size` must be at least 4.", call. = FALSE)
  }

  message("DISCLAIMER: This algorithm still in development by its authors.")

  matrix_profile <- matrix(Inf, matrix_profile_size, 1)
  left_matrix_profile <- right_matrix_profile <- matrix_profile
  profile_index <- matrix(-1, matrix_profile_size, 1)
  left_profile_index <- right_profile_index <- profile_index

  j <- 1
  order <- (exclusion_zone + 1):num_queries
  ssize <- min(s_size, length(order))
  order <- order[sample(seq_len(length(order)), size = ssize)]

  tictac <- Sys.time()

  if (verbose > 0) {
    pb <- utils::txtProgressBar(min = 1, max = ssize, style = 3, width = 80)
  }

  if (verbose > 0) {
    on.exit(close(pb))
  }
  if (verbose > 1) {
    on.exit(beep(sounds[[1]]), TRUE)
  }
  # anytime must return the result always
  on.exit(return(list(
    rmp = right_matrix_profile, rpi = right_profile_index,
    lmp = left_matrix_profile, lpi = left_profile_index,
    mp = Re(sqrt(as.complex(matrix_profile))), pi = profile_index
  )), TRUE)

  pre <- mass_pre(data, data_size, query, query_size, window_size = window_size)

  for (i in order) {
    j <- j + 1

    distance_profile <- diagonal_dist(
      data, i, data_size, window_size, num_queries, pre$data_mean,
      pre$data_sd
    )

    # distance_profile <- Re(sqrt(distance_profile))
    # distance_profile <- sqrt(abs(distance_profile))

    pos1 <- i:matrix_profile_size
    pos2 <- 1:(matrix_profile_size - i + 1)

    ind <- (matrix_profile[pos1] > distance_profile)
    profile_index[pos1[ind]] <- pos2[ind]
    matrix_profile[pos1[ind]] <- distance_profile[ind]
    ind <- (matrix_profile[pos2] > distance_profile)
    profile_index[pos2[ind]] <- pos1[ind]
    matrix_profile[pos2[ind]] <- distance_profile[ind]

    # matrix_profile[isSkip] <- Inf
    # profile_index[isSkip] <- 0

    # # anytime version
    # # left matrix_profile
    # ind <- (distance_profile[i:matrix_profile_size] < left_matrix_profile[i:matrix_profile_size])
    # ind <- c(rep(FALSE, (i - 1)), ind) # pad left
    # left_matrix_profile[ind] <- distance_profile[ind]
    # left_profile_index[which(ind)] <- i

    # # right matrix_profile
    # ind <- (distance_profile[1:i] < right_matrix_profile[1:i])
    # ind <- c(ind, rep(FALSE, matrix_profile_size - i)) # pad right
    # right_matrix_profile[ind] <- distance_profile[ind]
    # right_profile_index[which(ind)] <- i

    if (verbose > 0) {
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

diagonal_dist <- function(data, idx, data_size, window_size, mp_size, data_mean, data_sd) {
  data <- as.matrix(data)
  x_term <- matrix(1, mp_size - idx + 1, 1) *
    (t(data[idx:(idx + window_size - 1), 1, drop = FALSE]) %*% data[1:window_size, 1, drop = FALSE])[1]
  m_term <- data[idx:(mp_size - 1)] *
    data[1:(mp_size - idx)]
  a_term <- data[(idx + window_size):length(data)] *
    data[(window_size + 1):(data_size - idx + 1)]
  if (mp_size != idx) {
    x_term[2:length(x_term)] <- x_term[2:length(x_term)] - cumsum(m_term) + cumsum(a_term)
  }

  distance_profile <- (x_term - window_size * data_mean[idx:length(data_mean)] * data_mean[1:(mp_size - idx + 1)]) /
    (window_size * data_sd[idx:length(data_sd)] * data_sd[1:(mp_size - idx + 1)])
  distance_profile <- 2 * window_size * (1 - distance_profile)

  return(distance_profile)
}
