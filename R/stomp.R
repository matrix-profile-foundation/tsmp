#' Univariate STOMP algorithm
#'
#' Computes the Matrix Profile and Profile Index for Univariate Time Series.
#'
#' @details
#' The Matrix Profile, has the potential to revolutionize time series data mining because of its
#' generality, versatility, simplicity and scalability. In particular it has implications for time
#' series motif discovery, time series joins, shapelet discovery (classification), density
#' estimation, semantic segmentation, visualization, rule discovery, clustering etc. `verbose`
#' changes how much information is printed by this function; `0` means nothing, `1` means text, `2`
#' means text and sound. `exclusion.zone` is used to avoid  trivial matches; if a query data is
#' provided (join similarity), this parameter is ignored.
#'
#' @param ... a `matrix` or a `vector`. If a second time series is supplied it will be a join matrix
#'   profile.
#' @param window.size an `int`. Size of the sliding window.
#' @param exclusion.zone a `numeric`. Size of the exclusion zone, based on window size (default is
#'   `1/2`). See details.
#' @param verbose an `int`. See details. (Default is `2`).
#'
#' @return Returns the matrix profile `mp` and profile index `pi`. It also returns the left and
#'   right matrix profile `lmp`, `rmp` and profile index `lpi`, `rpi` that may be used to detect
#'   Time Series Chains (Yan Zhu 2018).
#' @export
#'
#' @family matrix profile computations
#' @seealso [stamp()], [stamp.par()]; [mstomp()], [mstomp.par()] for multivariate analysis.
#' @references * Zhu Y, Zimmerman Z, Senobari NS, Yeh CM, Funning G. Matrix Profile II : Exploiting
#'   a Novel Algorithm and GPUs to Break the One Hundred Million Barrier for Time Series Motifs and
#'   Joins. Icdm. 2016 Jan 22;54(1):739–48.
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#' mp <- stomp(toy_data$data[1:200,1], window.size = 30, verbose = 0)
#' \dontrun{
#' ref.data <- toy_data$data[,1]
#' query.data <- toy_data$data[,2]
#' # self similarity
#' mp <- stomp(ref.data, window.size = 30)
#' # join similarity
#' mp <- stomp(ref.data, query.data, window.size = 30)
#' }
stomp <- function(..., window_size, exclusion_zone = 1 / 2, verbose = 2) {
  args <- list(...)
  data <- args[[1]]
  if (length(args) > 1) {
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

  if (window_size > query_size / 2) {
    stop("Error: Time series is too short relative to desired window size.", call. = FALSE)
  }
  if (window_size < 4) {
    stop("Error: `window_size` must be at least 4.", call. = FALSE)
  }

  if (verbose > 0) {
    pb <- utils::txtProgressBar(min = 0, max = num_queries, style = 3, width = 80)
    on.exit(close(pb))
  }
  if (verbose > 1) {
    on.exit(beep(sounds[[1]]), TRUE)
  }

  data_fft <- matrix(0, (window_size + data_size), 1)
  data_mean <- matrix(0, matrix_profile_size, 1)
  data_sd <- matrix(0, matrix_profile_size, 1)
  first_product <- matrix(0, num_queries, 1)

  # forward
  nnpre <- mass_pre(data, data_size, query, query_size, window_size = window_size)
  data_fft <- nnpre$data_fft
  data_mean <- nnpre$data_mean
  data_sd <- nnpre$data_sd
  query_mean <- nnpre$query_mean
  query_sd <- nnpre$query_sd

  # reverse
  # This is needed to handle with the join similarity.
  rnnpre <- mass_pre(query, query_size, data, data_size, window_size = window_size)
  rdata_fft <- rnnpre$data_fft
  rdata_mean <- rnnpre$data_mean
  rdata_sd <- rnnpre$data_sd
  rquery_mean <- rnnpre$query_mean
  rquery_sd <- rnnpre$query_sd

  rnn <- mass(
    rdata_fft, data[1:window_size], query_size, window_size, rdata_mean,
    rdata_sd, rquery_mean[1], rquery_sd[1]
  )

  first_product[, 1] <- rnn$last_product

  tictac <- Sys.time()

  matrix_profile <- matrix(Inf, matrix_profile_size, 1)
  profile_index <- matrix(-1, matrix_profile_size, 1)
  left_matrix_profile <- matrix(Inf, matrix_profile_size, 1)
  left_profile_index <- matrix(-1, matrix_profile_size, 1)
  right_matrix_profile <- matrix(Inf, matrix_profile_size, 1)
  right_profile_index <- matrix(-1, matrix_profile_size, 1)
  distance_profile <- matrix(0, matrix_profile_size, 1)
  last_product <- matrix(0, matrix_profile_size, 1)
  drop_value <- matrix(0, 1, 1)

  for (i in 1:num_queries) {
    # compute the distance profile
    query_window <- as.matrix(query[i:(i + window_size - 1), 1])

    if (i == 1) {
      nn <- mass(
        data_fft, query_window, data_size, window_size, data_mean, data_sd,
        query_mean[i], query_sd[i]
      )
      distance_profile[, 1] <- nn$distance_profile
      last_product[, 1] <- nn$last_product
    } else {
      last_product[2:(data_size - window_size + 1), 1] <- last_product[1:(data_size - window_size), 1] -
        data[1:(data_size - window_size), 1] * drop_value +
        data[(window_size + 1):data_size, 1] * query_window[window_size, 1]

      last_product[1, 1] <- first_product[i, 1]
      distance_profile <- 2 * (window_size - (last_product - window_size * data_mean * query_mean[i]) /
        (data_sd * query_sd[i]))
    }

    distance_profile <- Re(sqrt(distance_profile))
    drop_value <- query_window[1, 1]

    # apply exclusion zone
    if (exclusion_zone > 0) {
      exc_st <- max(1, i - exclusion_zone)
      exc_ed <- min(matrix_profile_size, i + exclusion_zone)
      distance_profile[exc_st:exc_ed, 1] <- Inf
      distance_profile[data_sd < vars()$eps] <- Inf
      if (skip_location[i] || any(query_sd[i] < vars()$eps)) {
        distance_profile[] <- Inf
      }
    }

    # left matrix_profile
    ind <- (distance_profile[i:matrix_profile_size] < left_matrix_profile[i:matrix_profile_size])
    ind <- c(rep(FALSE, (i - 1)), ind) # pad left
    left_matrix_profile[ind] <- distance_profile[ind]
    left_profile_index[which(ind)] <- i

    # right matrix_profile
    ind <- (distance_profile[1:i] < right_matrix_profile[1:i])
    ind <- c(ind, rep(FALSE, matrix_profile_size - i)) # pad right
    right_matrix_profile[ind] <- distance_profile[ind]
    right_profile_index[which(ind)] <- i

    ind <- (distance_profile < matrix_profile)
    matrix_profile[ind] <- distance_profile[ind]
    profile_index[which(ind)] <- i

    if (verbose > 0) {
      utils::setTxtProgressBar(pb, i)
    }
  }

  tictac <- Sys.time() - tictac

  if (verbose > 0) {
    message(sprintf("\nFinished in %.2f %s", tictac, units(tictac)))
  }

  return(list(
    rmp = right_matrix_profile, rpi = right_profile_index,
    lmp = left_matrix_profile, lpi = left_profile_index,
    mp = matrix_profile, pi = profile_index
  ))
}
