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
#' adds the progress bar, `3` adds the finish sound. `exclusion_zone` is used to avoid  trivial
#' matches; if a query data is provided (join similarity), this parameter is ignored.
#'
#' @param \dots a `matrix` or a `vector`. If a second time series is supplied it will be a join matrix
#'   profile.
#' @param window_size an `int`. Size of the sliding window.
#' @param exclusion_zone a `numeric`. Size of the exclusion zone, based on window size (default is
#'   `1/2`). See details.
#' @param verbose an `int`. See details. (Default is `2`).
#'
#' @return Returns a `MatrixProfile` object, a `list` with the matrix profile `mp`, profile index `pi`
#'   left and right matrix profile `lmp`, `rmp` and profile index `lpi`, `rpi`, window size `w` and
#'   exclusion zone `ez`.
#'
#' @export
#'
#' @family matrix profile computations
#'
#' @describeIn stomp Single thread version.
#'
#' @references * Zhu Y, Zimmerman Z, Senobari NS, Yeh CM, Funning G. Matrix Profile II : Exploiting
#'   a Novel Algorithm and GPUs to Break the One Hundred Million Barrier for Time Series Motifs and
#'   Joins. Icdm. 2016 Jan 22;54(1):739-48.
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#' mp <- stomp(mp_toy_data$data[1:200, 1], window_size = 30, verbose = 0)
#' \dontrun{
#' #' # using threads
#' mp <- stomp_par(mp_toy_data$data[1:400, 1], window_size = 30, verbose = 0)
#'
#' ref_data <- mp_toy_data$data[, 1]
#' query_data <- mp_toy_data$data[, 2]
#' # self similarity
#' mp <- stomp(ref_data, window_size = 30)
#' # join similarity
#' mp2 <- stomp(ref_data, query_data, window_size = 30)
#' }
stomp <- function(..., window_size, exclusion_zone = getOption("tsmp.exclusion_zone", 1 / 2),
                  verbose = getOption("tsmp.verbose", 2)) {
  argv <- list(...)
  argc <- length(argv)
  data <- argv[[1]]
  if (argc > 1 && !is.null(argv[[2]])) {
    query <- argv[[2]]
    exclusion_zone <- 0 # don't use exclusion zone for joins
    join <- TRUE
  } else {
    query <- data
    join <- FALSE
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
    stop("Unknown type of data. Must be: a column matrix or a vector.", call. = FALSE)
  }

  if (is.vector(query)) {
    query <- as.matrix(query)
  } else if (is.matrix(query)) {
    if (ncol(query) > nrow(query)) {
      query <- t(query)
    }
  } else {
    stop("Unknown type of query. Must be: a column matrix or a vector.", call. = FALSE)
  }

  ez <- exclusion_zone # store original
  exclusion_zone <- round(window_size * exclusion_zone + vars()$eps)
  data_size <- nrow(data)
  query_size <- nrow(query)
  matrix_profile_size <- data_size - window_size + 1
  num_queries <- query_size - window_size + 1

  if (query_size > data_size) {
    stop("Query must be smaller or the same size as reference data.")
  }
  if (window_size > ceiling(query_size / 2)) {
    stop("Time series is too short relative to desired window size.", call. = FALSE)
  }
  if (window_size < 4) {
    stop("`window_size` must be at least 4.", call. = FALSE)
  }

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

  if (verbose > 1) {
    pb <- progress::progress_bar$new(
      format = "STOMP [:bar] :percent at :tick_rate it/s, elapsed: :elapsed, eta: :eta",
      clear = FALSE, total = num_queries, width = 80
    )
  }

  if (verbose > 2) {
    on.exit(beep(sounds[[1]]), TRUE)
  }

  first_product <- matrix(0, num_queries, 1)

  # forward
  nn <- dist_profile(data, query, window_size = window_size)
  # reverse
  # This is needed to handle with the join similarity.
  rnn <- dist_profile(query, data, window_size = window_size)

  first_product[, 1] <- rnn$last_product

  tictac <- Sys.time()

  matrix_profile <- matrix(Inf, matrix_profile_size, 1)
  profile_index <- matrix(-Inf, matrix_profile_size, 1)
  if (join) {
    # no RMP and LMP for joins
    left_matrix_profile <- right_matrix_profile <- NULL
    left_profile_index <- right_profile_index <- NULL
  } else {
    left_matrix_profile <- right_matrix_profile <- matrix_profile
    left_profile_index <- right_profile_index <- profile_index
  }
  distance_profile <- matrix(0, matrix_profile_size, 1)
  last_product <- matrix(0, matrix_profile_size, 1)
  drop_value <- matrix(0, 1, 1)

  for (i in 1:num_queries) {
    # compute the distance profile
    query_window <- as.matrix(query[i:(i + window_size - 1), 1])

    if (i == 1) {
      distance_profile[, 1] <- nn$distance_profile
      last_product[, 1] <- nn$last_product
    } else {
      last_product[2:(data_size - window_size + 1), 1] <- last_product[1:(data_size - window_size), 1] -
        data[1:(data_size - window_size), 1] * drop_value +
        data[(window_size + 1):data_size, 1] * query_window[window_size, 1]

      last_product[1, 1] <- first_product[i, 1]
      distance_profile <- 2 * (window_size - (last_product - window_size * nn$par$data_mean * nn$par$query_mean[i]) /
        (nn$par$data_sd * nn$par$query_sd[i]))
    }

    distance_profile[distance_profile < 0] <- 0
    distance_profile <- sqrt(distance_profile)
    drop_value <- query_window[1, 1]

    # apply exclusion zone
    if (exclusion_zone > 0) {
      exc_st <- max(1, i - exclusion_zone)
      exc_ed <- min(matrix_profile_size, i + exclusion_zone)
      distance_profile[exc_st:exc_ed, 1] <- Inf
    }

    distance_profile[nn$par$data_sd < vars()$eps] <- Inf
    if (skip_location[i] || any(nn$par$query_sd[i] < vars()$eps)) {
      distance_profile[] <- Inf
    }
    distance_profile[skip_location] <- Inf

    if (!join) {
      # no RMP and LMP for joins
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
    }

    ind <- (distance_profile < matrix_profile)
    matrix_profile[ind] <- distance_profile[ind]
    profile_index[which(ind)] <- i

    if (verbose > 1) {
      pb$tick()
    }
  }

  tictac <- Sys.time() - tictac

  if (verbose > 0) {
    message(sprintf("Finished in %.2f %s", tictac, units(tictac)))
  }

  return({
    obj <- list(
      mp = matrix_profile, pi = profile_index,
      rmp = right_matrix_profile, rpi = right_profile_index,
      lmp = left_matrix_profile, lpi = left_profile_index,
      w = window_size,
      ez = ez
    )
    class(obj) <- "MatrixProfile"
    attr(obj, "join") <- join
    obj
  })
}
