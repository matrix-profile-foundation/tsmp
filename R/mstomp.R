#' Multivariate STOMP algorithm
#'
#' Computes the Matrix Profile and Profile Index for Multivariate Time Series.
#'
#' @details
#' The Matrix Profile, has the potential to revolutionize time series data mining because of its
#' generality, versatility, simplicity and scalability. In particular it has implications for time
#' series motif discovery, time series joins, shapelet discovery (classification), density
#' estimation, semantic segmentation, visualization, rule discovery, clustering etc. The MSTOMP
#' computes the Matrix Profile and Profile Index for Multivariate Time Series that is meaningful for
#' multidimensional MOTIF discovery. It uses the STOMP algorithm that is faster than STAMP but lacks
#' its anytime property.
#'
#' Although this functions handles Multivariate Time Series, it can also be used to handle
#' Univariate Time Series. `verbose` changes how much information is printed by this function; `0`
#' means nothing, `1` means text, `2` adds the progress bar, `3` adds the finish sound.
#'
#' @param data a `matrix` of `numeric`, where each column is a time series. Accepts `vector` (see
#'   details), `list` and `data.frame` too.
#' @param window_size an `int` with the size of the sliding window.
#' @param exclusion_zone a `numeric`. Size of the exclusion zone, based on window size (default is
#'   `1/2`).
#' @param verbose an `int`. See details. (Default is `2`).
#' @param must_dim an `int` or `vector` of which dimensions to forcibly include (default is `NULL`).
#' @param exc_dim an `int` or `vector` of which dimensions to exclude (default is `NULL`).
#'
#' @return Returns a `MultiMatrixProfile` object, a `list` with the matrix profile `mp`, profile index `pi`
#'   left and right matrix profile `lmp`, `rmp` and profile index `lpi`, `rpi`, window size `w`,
#'   number of dimensions `n_dim`, exclusion zone `ez`, must dimensions `must` and excluded dimensions `exc`.
#'
#'   If the input has only one dimension, returns the same as [stomp()].
#'
#' @export
#'
#' @describeIn mstomp Single thread version.
#'
#' @family matrix profile computations
#'
#' @references * Yeh CM, Kavantzas N, Keogh E. Matrix Profile VI : Meaningful Multidimensional Motif
#'   Discovery.
#' @references * Zhu Y, Imamura M, Nikovski D, Keogh E. Matrix Profile VII: Time Series Chains: A
#'   New Primitive for Time Series Data Mining. Knowl Inf Syst. 2018 Jun 2;1-27.
#' @references Website: <https://sites.google.com/view/mstamp/>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#' # using all dimensions
#' mp <- mstomp(mp_toy_data$data[1:150, ], 30, verbose = 0)
#' \dontrun{
#' #' # using threads
#' mp <- mstomp_par(mp_toy_data$data[1:150, ], 30, verbose = 0)
#'
#' # force using dimensions 1 and 2
#' mp <- mstomp(mp_toy_data$data[1:200, ], 30, must_dim = c(1, 2))
#' # exclude dimensions 2 and 3
#' mp2 <- mstomp(mp_toy_data$data[1:200, ], 30, exc_dim = c(2, 3))
#' }
#'
mstomp <- function(data, window_size, exclusion_zone = getOption("tsmp.exclusion_zone", 1 / 2),
                   verbose = getOption("tsmp.verbose", 2),
                   must_dim = NULL, exc_dim = NULL) {
  # get various length
  ez <- exclusion_zone # store original
  exclusion_zone <- round(window_size * exclusion_zone + vars()$eps)

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

  matrix_profile_size <- data_size - window_size + 1

  # check input
  if (window_size > data_size / 2) {
    stop("Time series is too short relative to desired window size.")
  }
  if (window_size < 4) {
    stop("`window_size` must be at least 4.")
  }
  if (any(must_dim > n_dim)) {
    stop("`must_dim` must be less then the total dimension.")
  }
  if (any(exc_dim > n_dim)) {
    stop("`exc_dim` must be less then the total dimension.")
  }
  if (length(intersect(must_dim, exc_dim)) > 0) {
    stop("The same dimension is presented in both the exclusion dimension and must have dimension.")
  }

  # check skip position
  n_exc <- length(exc_dim)
  n_must <- length(must_dim)
  mask_exc <- rep(FALSE, n_dim)
  mask_exc[exc_dim] <- TRUE
  skip_location <- rep(FALSE, matrix_profile_size)

  for (i in 1:matrix_profile_size) {
    if (any(is.na(data[i:(i + window_size - 1), !mask_exc])) || any(is.infinite(data[i:(i + window_size - 1), !mask_exc]))) {
      skip_location[i] <- TRUE
    }
  }

  data[is.na(data)] <- 0
  data[is.infinite(data)] <- 0

  if (verbose > 1) {
    pb <- progress::progress_bar$new(
      format = "mSTOMP [:bar] :percent at :tick_rate it/s, elapsed: :elapsed, eta: :eta",
      clear = FALSE, total = matrix_profile_size, width = 80
    )
  }
  if (verbose > 2) {
    on.exit(beep(sounds[[1]]), TRUE)
  }

  # initialization
  nn <- vector(mode = "list", length = 3)
  data_mean <- matrix(0, matrix_profile_size, n_dim)
  data_sd <- matrix(0, matrix_profile_size, n_dim)
  first_product <- matrix(0, matrix_profile_size, n_dim)


  for (i in 1:n_dim) {
    nn[[i]] <- dist_profile(data[, i], data[, i], window_size = window_size)
    first_product[, i] <- nn[[i]]$last_product
    data_mean[, i] <- nn[[i]]$par$data_mean
    data_sd[, i] <- nn[[i]]$par$data_sd
  }

  tictac <- Sys.time()
  # compute the matrix profile
  matrix_profile <- matrix(Inf, matrix_profile_size, n_dim)
  profile_index <- matrix(-Inf, matrix_profile_size, n_dim)
  left_matrix_profile <- matrix(Inf, matrix_profile_size, n_dim)
  left_profile_index <- matrix(-Inf, matrix_profile_size, n_dim)
  right_matrix_profile <- matrix(Inf, matrix_profile_size, n_dim)
  right_profile_index <- matrix(-Inf, matrix_profile_size, n_dim)
  distance_profile <- matrix(0, matrix_profile_size, n_dim)
  last_product <- matrix(0, matrix_profile_size, n_dim)
  drop_value <- matrix(0, 1, n_dim)

  for (i in 1:matrix_profile_size) {
    # compute the distance profile
    if (verbose > 1) {
      pb$tick()
    }

    query_window <- as.matrix(data[i:(i + window_size - 1), ])

    if (i == 1) {
      for (j in 1:n_dim) {
        distance_profile[, j] <- nn[[j]]$distance_profile
        last_product[, j] <- nn[[j]]$last_product
      }
    } else {
      rep_drop_value <- kronecker(matrix(1, matrix_profile_size - 1, 1), t(drop_value))
      rep_query <- kronecker(matrix(1, matrix_profile_size - 1, 1), t(query_window[window_size, ]))

      last_product[2:(data_size - window_size + 1), ] <- last_product[1:(data_size - window_size), ] -
        data[1:(data_size - window_size), ] * rep_drop_value +
        data[(window_size + 1):data_size, ] * rep_query


      last_product[1, ] <- first_product[i, ]

      distance_profile <- 2 * (window_size - (last_product - window_size * data_mean * kronecker(matrix(1, matrix_profile_size, 1), t(data_mean[i, ]))) /
        (data_sd * kronecker(matrix(1, matrix_profile_size, 1), t(data_sd[i, ]))))
    }

    drop_value <- query_window[1, ]

    # apply exclusion zone
    exc_st <- max(1, i - exclusion_zone)
    exc_ed <- min(matrix_profile_size, i + exclusion_zone)
    distance_profile[exc_st:exc_ed, ] <- Inf
    distance_profile[data_sd < vars()$eps] <- Inf
    if (skip_location[i] || any(data_sd[i, !mask_exc] < vars()$eps)) {
      distance_profile[] <- Inf
    }
    distance_profile[skip_location, ] <- Inf

    # apply dimension "must have" and "exclusion"
    distance_profile[, exc_dim] <- Inf

    if (n_must > 0) {
      mask_must <- rep(FALSE, n_dim)
      mask_must[must_dim] <- TRUE
      dist_pro_must <- distance_profile[, mask_must]
      distance_profile[, mask_must] <- -Inf
    }

    if (n_dim > 1) {
      dist_pro_sort <- t(apply(distance_profile, 1, sort))
    } # sort by row, put all -Inf to the first column
    else {
      dist_pro_sort <- distance_profile
    }

    if (n_must > 0) {
      dist_pro_sort[, 1:n_must] <- dist_pro_must
    }

    # figure out and store the nearest neighbor
    dist_pro_cum <- rep(0, matrix_profile_size)
    dist_pro_merg <- rep(0, matrix_profile_size)

    for (j in (max(1, n_must):(n_dim - n_exc))) {
      dist_pro_cum <- dist_pro_cum + dist_pro_sort[, j]
      dist_pro_merg[] <- dist_pro_cum / j

      # left matrix_profile
      if (i > (exclusion_zone + 1)) {
        min_idx <- which.min(dist_pro_merg[1:(i - exclusion_zone)])
        min_val <- dist_pro_merg[min_idx]
        left_matrix_profile[i, j] <- min_val
        left_profile_index[i, j] <- min_idx
      }

      # right matrix_profile
      if (i < (matrix_profile_size - exclusion_zone)) {
        min_idx <- which.min(dist_pro_merg[(i + exclusion_zone):matrix_profile_size]) + i + exclusion_zone - 1
        min_val <- dist_pro_merg[min_idx]
        right_matrix_profile[i, j] <- min_val
        right_profile_index[i, j] <- min_idx
      }

      # normal matrix_profile
      min_idx <- which.min(dist_pro_merg)
      min_val <- dist_pro_merg[min_idx]
      matrix_profile[i, j] <- min_val
      profile_index[i, j] <- min_idx
    }
  }

  matrix_profile <- sqrt(matrix_profile)
  right_matrix_profile <- sqrt(right_matrix_profile)
  left_matrix_profile <- sqrt(left_matrix_profile)

  # remove bad k setting in the returned matrix
  if (n_must > 1) {
    matrix_profile[, 1:(n_must - 1)] <- NA
    right_matrix_profile[, 1:(n_must - 1)] <- NA
    left_matrix_profile[, 1:(n_must - 1)] <- NA
  }
  if (n_exc > 0) {
    matrix_profile[, (n_dim - n_exc + 1):n_dim] <- NA
    right_matrix_profile[, (n_dim - n_exc + 1):n_dim] <- NA
    left_matrix_profile[, (n_dim - n_exc + 1):n_dim] <- NA
  }
  if (n_must > 1) {
    profile_index[, 1:(n_must - 1)] <- NA
    right_profile_index[, 1:(n_must - 1)] <- NA
    left_profile_index[, 1:(n_must - 1)] <- NA
  }
  if (n_exc > 0) {
    profile_index[, (n_dim - n_exc + 1):n_dim] <- NA
    right_profile_index[, (n_dim - n_exc + 1):n_dim] <- NA
    left_profile_index[, (n_dim - n_exc + 1):n_dim] <- NA
  }

  tictac <- Sys.time() - tictac

  if (verbose > 0) {
    message(sprintf("Finished in %.2f %s", tictac, units(tictac)))
  }

  if (n_dim > 1) {
    obj <- list(
      mp = matrix_profile, pi = profile_index,
      rmp = right_matrix_profile, rpi = right_profile_index,
      lmp = left_matrix_profile, lpi = left_profile_index,
      w = window_size,
      ez = ez,
      n_dim = n_dim,
      must = must_dim,
      exc = exc_dim
    )
    class(obj) <- "MultiMatrixProfile"
    attr(obj, "join") <- FALSE
  } else {
    obj <- list(
      mp = matrix_profile, pi = profile_index,
      rmp = right_matrix_profile, rpi = right_profile_index,
      lmp = left_matrix_profile, lpi = left_profile_index,
      w = window_size,
      ez = ez
    )
    class(obj) <- "MatrixProfile"
    attr(obj, "join") <- FALSE
  }

  return(obj)
}
