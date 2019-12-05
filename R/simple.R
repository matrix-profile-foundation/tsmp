#' Compute the join similarity for Sound data
#'
#' @details
#' `verbose` changes how much information is printed by this function; `0` means nothing, `1` means
#' text, `2` adds the progress bar, `3` adds the finish sound.
#'
#' @param \dots a `matrix` of `numeric`, where each column is a time series. Accepts `list` and
#'   `data.frame` too. If a second time series is supplied it will be a join matrix profile.
#' @param window_size an `int` with the size of the sliding window.
#' @param exclusion_zone a `numeric`. Size of the exclusion zone, based on window size (default is
#'   `1/2`).
#' @param verbose an `int`. See details. (Default is `2`).
#'
#' @return Returns a `SimpleMatrixProfile` object, a `list` with the matrix profile `mp`, profile index `pi`,
#' number of dimensions `n_dim`, window size `w` and exclusion zone `ez`.
#'
#' @export
#' @references * Silva D, Yeh C, Batista G, Keogh E. Simple: Assessing Music Similarity Using
#'   Subsequences Joins. Proc 17th ISMIR Conf. 2016;23-30.
#' @references * Silva DF, Yeh C-CM, Zhu Y, Batista G, Keogh E. Fast Similarity Matrix Profile for
#'   Music Analysis and Exploration. IEEE Trans Multimed. 2018;14(8):1-1.
#' @references Website: <https://sites.google.com/view/simple-fast>
#' @references Website: <https://sites.google.com/site/ismir2016simple/home>
#'
#' @examples
#' w <- 30
#' data <- mp_toy_data$data # 3 dimensions matrix
#' result <- simple_fast(data, window_size = w, verbose = 0)
simple_fast <- function(..., window_size, exclusion_zone = getOption("tsmp.exclusion_zone", 1 / 2),
                        verbose = getOption("tsmp.verbose", 2)) {
  if (!is.numeric(window_size) || length(window_size) > 1) {
    stop("Unknown type of `window_size`. Must be an `int` or `numeric`")
  }
  argv <- list(...)
  argc <- length(argv)
  data <- argv[[1]]
  join <- FALSE
  if (argc > 1 && !is.null(argv[[2]])) {
    query <- argv[[2]]
    if (!isTRUE(all.equal(data, query))) {
      exclusion_zone <- 0
      join <- TRUE
    } # don't use exclusion zone for joins
  } else {
    query <- data
  }

  # transform data list into matrix
  if (is.matrix(data) || is.data.frame(data)) {
    if (is.data.frame(data)) {
      data <- as.matrix(data)
    } # just to be uniform
    if (ncol(data) > nrow(data)) {
      data <- t(data)
    }
    data_size <- nrow(data)
    n_dim <- ncol(data)
  } else if (is.list(data)) {
    data_size <- length(data[[1]])
    n_dim <- length(data)

    for (i in 1:n_dim) {
      len <- length(data[[i]])
      # Fix TS size with NaN
      if (len < data_size) {
        data[[i]] <- c(data[[i]], rep(NA, data_size - len))
      }
    }
    # transform data into matrix (each column is a TS)
    data <- sapply(data, cbind)
  } else if (is.vector(data)) {
    data_size <- length(data)
    n_dim <- 1
    # transform data into 1-col matrix
    data <- as.matrix(data) # just to be uniform
  } else {
    stop("Unknown type of data. Must be: matrix, data.frame, vector or list.")
  }

  # transform query list into matrix
  if (is.matrix(query) || is.data.frame(query)) {
    if (is.data.frame(query)) {
      query <- as.matrix(query)
    } # just to be uniform
    if (ncol(query) > nrow(query)) {
      query <- t(query)
    }
    query_size <- nrow(query)
    q_dim <- ncol(query)
  } else if (is.list(query)) {
    query_size <- length(query[[1]])
    q_dim <- length(query)

    for (i in 1:q_dim) {
      len <- length(query[[i]])
      # Fix TS size with NaN
      if (len < query_size) {
        query[[i]] <- c(query[[i]], rep(NA, query_size - len))
      }
    }
    # transform query into matrix (each column is a TS)
    query <- sapply(query, cbind)
  } else if (is.vector(query)) {
    query_size <- length(query)
    q_dim <- 1
    # transform query into 1-col matrix
    query <- as.matrix(query) # just to be uniform
  } else {
    stop("Unknown type of query. Must be: matrix, data.frame, vector or list.")
  }

  # check input
  if (q_dim != n_dim) {
    stop("Data and query dimensions must be the same.")
  }
  if (window_size > data_size / 2) {
    stop("Reference Time series is too short relative to desired window size.")
  }
  if (window_size > query_size / 2) {
    stop("Query Time series is too short relative to desired window size.")
  }
  if (window_size < 4) {
    stop("`window_size` must be at least 4.")
  }

  ez <- exclusion_zone # store original
  exclusion_zone <- round(window_size * exclusion_zone + vars()$eps)

  # initialization
  matrix_profile_size <- data_size - window_size + 1
  matrix_profile <- rep(Inf, matrix_profile_size)
  profile_index <- rep(-Inf, matrix_profile_size)

  if (verbose > 1) {
    pb <- progress::progress_bar$new(
      format = "SiMPle [:bar] :percent at :tick_rate it/s, elapsed: :elapsed, eta: :eta",
      clear = FALSE, total = matrix_profile_size, width = 80
    )
  }
  if (verbose > 2) {
    on.exit(beep(sounds[[1]]), TRUE)
  }

  # for the first dot-product for both data and query
  pre_data <- mass_simple_pre(data, data_size, window_size = window_size)
  data_fft <- pre_data$data_fft
  data_sumx2 <- pre_data$sumx2
  query_window <- query[1:window_size, ]
  res_data <- mass_simple(data_fft, query_window, data_size, window_size, data_sumx2)

  if (verbose > 1) {
    pb$tick()
  }

  pre_query <- mass_simple_pre(query, query_size, window_size = window_size)
  query_fft <- pre_query$data_fft
  query_sumx2 <- pre_query$sumx2
  data_window <- data[1:window_size, ]
  res_query <- mass_simple(query_fft, data_window, query_size, window_size, query_sumx2)

  distance_profile <- res_query$distance_profile
  last_product <- res_query$last_product
  first_product <- res_data$last_product
  query_sumy2 <- res_query$sumy2
  dropval <- data_window[1, ] # dropval is the first element of refdata window

  # no ez if join
  distance_profile[1:exclusion_zone] <- Inf


  ind <- which.min(distance_profile)
  profile_index[1] <- ind
  matrix_profile[1] <- distance_profile[ind]

  tictac <- Sys.time()

  # compute the remainder of the matrix profile
  for (i in 2:matrix_profile_size) {

    # compute the distance profile
    if (verbose > 1) {
      pb$tick()
    }

    data_window <- data[i:(i + window_size - 1), ]

    query_sumy2 <- query_sumy2 - dropval^2 + data_window[window_size, ]^2

    for (j in 1:n_dim) {
      last_product[2:(data_size - window_size + 1), j] <- last_product[1:(data_size - window_size), j] -
        query[1:(data_size - window_size), j] * dropval[j] +
        query[(window_size + 1):data_size, j] * data_window[window_size, j]
    }

    last_product[1, ] <- first_product[i, ]
    dropval <- data_window[1, ]

    distance_profile <- matrix(0, nrow(query_sumx2), 1)

    for (j in 1:n_dim) {
      distance_profile <- distance_profile + query_sumx2[, j] - 2 * last_product[, j] + query_sumy2[j]
    }

    if (exclusion_zone > 0) {
      exc_st <- max(1, i - exclusion_zone)
      exc_ed <- min(matrix_profile_size, i + exclusion_zone)
      distance_profile[exc_st:exc_ed] <- Inf
    }

    ind <- which.min(distance_profile)
    profile_index[i] <- ind
    matrix_profile[i] <- distance_profile[ind]
  }

  tictac <- Sys.time() - tictac

  if (verbose > 0) {
    message(sprintf("Finished in %.2f %s", tictac, units(tictac)))
  }

  obj <- list(
    mp = as.matrix(matrix_profile), pi = as.matrix(profile_index),
    rmp = NULL, rpi = NULL,
    lmp = NULL, lpi = NULL,
    n_dim = n_dim,
    w = window_size,
    ez = ez
  )
  class(obj) <- "SimpleMatrixProfile"
  attr(obj, "join") <- join
  return(obj)
}

#' Precomputes several values used on MASS
#'
#' The difference of this function to [mass_pre()] is that this does not normalize data. Specific for this domain.
#'
#' @param data a `matrix` of `numeric`. Reference Time Series.
#' @param data_size an `int`. Reference Time Series size.
#' @param window_size an `int`. Sliding window size.
#'
#' @return Returns `data_fft` and `sumx2`.
#' @keywords internal
#' @noRd
#'
#' @references * Abdullah Mueen, Yan Zhu, Michael Yeh, Kaveh Kamgar, Krishnamurthy Viswanathan, Chetan Kumar Gupta and Eamonn Keogh (2015), The Fastest Similarity Search Algorithm for Time Series Subsequences under Euclidean Distance.
#' @references Website: <https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html>


mass_simple_pre <- function(data, data_size, window_size) {
  if (nrow(data) < ncol(data)) {
    data <- t(data)
  }

  n_dim <- ncol(data)

  data <- rbind(data, matrix(0, data_size, n_dim))

  data_fft <- apply(data, 2, stats::fft)
  cum_sumx2 <- apply(data^2, 2, cumsum)

  sumx2 <- cum_sumx2[window_size:data_size, ] - rbind(rep(0, n_dim), cum_sumx2[1:(data_size - window_size), ])

  return(list(data_fft = data_fft, sumx2 = sumx2))
}

#' Calculates the distance profile using MASS algorithm
#'
#' Mueen's Algorithm for Similarity Search is The Fastest Similarity Search Algorithm for Time Series Subsequences under Euclidean Distance and Correlation Coefficient.
#' The difference of this function to [mass()] is that this does not normalize data. Specific for this domain.
#'
#' @param data_fft precomputed data product.
#' @param query_window a `matrix` of `numeric`. Query window.
#' @param data_size an `int`. The length of the reference data.
#' @param window_size an `int`. Sliding window size.
#' @param sumx2 precomputed sum of squares
#'
#' @return Returns the `distance_profile` for the given query and the `last_product` for STOMP algorithm and `sumy2`.
#' @keywords internal
#' @noRd
#'
#' @references * Abdullah Mueen, Yan Zhu, Michael Yeh, Kaveh Kamgar, Krishnamurthy Viswanathan, Chetan Kumar Gupta and Eamonn Keogh (2015), The Fastest Similarity Search Algorithm for Time Series Subsequences under Euclidean Distance
#' @references Website: <https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html>
#'


mass_simple <- function(data_fft, query_window, data_size, window_size, sumx2) {
  if (nrow(data_fft) < ncol(data_fft)) {
    data_fft <- t(data_fft)
  }

  n_dim <- ncol(data_fft)

  if (ncol(query_window) != n_dim) {
    query_window <- t(query_window)
  }

  # pre-process query for fft
  query_window <- apply(query_window, 2, rev)
  query_window <- rbind(query_window, matrix(0, 2 * data_size - window_size, n_dim))

  query_fft <- apply(query_window, 2, stats::fft)
  # compute the product
  prod <- data_fft * query_fft
  z <- apply(prod, 2, function(x) {
    stats::fft(x, inverse = TRUE) / length(x)
  })

  sumy2 <- apply(query_window^2, 2, sum)

  last_product <- Re(z[window_size:data_size, ])

  distance_profile <- matrix(0, nrow(sumx2), 1)

  for (i in 1:n_dim) {
    distance_profile <- distance_profile + sumx2[, i] - 2 * last_product[, i] + sumy2[i]
  }

  return(list(distance_profile = distance_profile, last_product = last_product, sumy2 = sumy2))
}
