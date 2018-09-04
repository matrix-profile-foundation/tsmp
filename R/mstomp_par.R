#' Multivariate STOMP algorithm Parallel version
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
#' means nothing, `1` means text, `2` means text and sound.
#'
#' @param data a `matrix` of `numeric`, where each column is a time series. Accepts `vector` (see
#'   details), `list` and `data.frame` too.
#' @param window.size an `int`. Size of the sliding window.
#' @param must.dim an `int` or `vector` of which dimensions to forcibly include (default is `NULL`).
#' @param exc.dim an `int` or `vector` of which dimensions to exclude (default is `NULL`).
#' @param exclusion.zone a `numeric`. Size of the exclusion zone, based on window size (default is
#'   `1/2`). See details.
#' @param verbose an `int`. See details. (Default is `2`).
#' @param n.workers an `int`. Number of workers for parallel. (Default is `2`).
#'
#' @return Returns the matrix profile `mp` and profile index `pi`. It also returns the left and
#'   right matrix profile `lmp`, `rmp` and profile index `lpi`, `rpi` that may be used to detect
#'   Time Series Chains (Yan Zhu 2018).
#' @export
#'
#' @family matrix profile computations
#' @seealso [stamp()], [stamp.par()], [mstomp()]
#' @references * Yeh CM, Kavantzas N, Keogh E. Matrix Profile VI : Meaningful Multidimensional Motif
#'   Discovery.
#' @references * Zhu Y, Imamura M, Nikovski D, Keogh E. Matrix Profile VII: Time Series Chains: A
#'   New Primitive for Time Series Data Mining. Knowl Inf Syst. 2018 Jun 2;1–27.
#' @references Website: <https://sites.google.com/view/mstamp/>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#' # using all dimensions
#' mp <- mstomp.par(toy_data$data[1:100,], 30, verbose = 0)
#' @import doSNOW foreach parallel

mstomp_par <- function(data, window_size, must_dim = NULL, exc_dim = NULL, exclusion_zone = 1 / 2, verbose = 2, n_workers = 2) {
  ## get various length
  exclusion_zone <- floor(window_size * exclusion_zone)

  ## transform data list into matrix
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
    stop("Error: Unknown type of data. Must be: matrix, data_frame, vector or list.", call. = FALSE)
  }

  matrix_profile_size <- data_size - window_size + 1

  ## check input
  if (window_size > data_size / 2) {
    stop("Error: Time series is too short relative to desired window size.", call. = FALSE)
  }
  if (window_size < 4) {
    stop("Error: `window_size` must be at least 4.", call. = FALSE)
  }
  if (any(must_dim > n_dim)) {
    stop("Error: `must_dim` must be less then the total dimension.", call. = FALSE)
  }
  if (any(exc_dim > n_dim)) {
    stop("Error: `exc_dim` must be less then the total dimension.", call. = FALSE)
  }
  if (length(intersect(must_dim, exc_dim)) > 0) {
    stop("Error: The same dimension is presented in both the exclusion dimension and must have dimension.", call. = FALSE)
  }

  ## check skip position
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

  ## initialization
  data_fft <- matrix(0, (window_size + data_size), n_dim)
  data_mean <- matrix(0, matrix_profile_size, n_dim)
  data_sd <- matrix(0, matrix_profile_size, n_dim)
  first_product <- matrix(0, matrix_profile_size, n_dim)

  for (i in 1:n_dim) {
    nnpre <- mass_pre(data[, i], data_size, window_size = window_size)
    data_fft[, i] <- nnpre$data_fft
    data_mean[, i] <- nnpre$data_mean
    data_sd[, i] <- nnpre$data_sd

    mstomp <- mass(data_fft[, i], data[1:window_size, i], data_size, window_size, data_mean[, i], data_sd[, i], data_mean[1, i], data_sd[1, i])
    first_product[, i] <- mstomp$last_product
  }

  cores <- min(max(2, n_workers), parallel::detectCores())

  # SNOW package
  if (verbose > 0) {
    progress <- function(n) utils::setTxtProgressBar(pb, n)
  }
  else {
    progress <- function(n) return(invisible(TRUE))
  }
  opts <- list(progress = progress)

  cl <- parallel::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  on.exit(parallel::stopCluster(cl))
  if (verbose > 0) {
    on.exit(close(pb), TRUE)
  }
  if (verbose > 1) {
    on.exit(beep(sounds[[1]]), TRUE)
  }

  ## initialize variable
  per_work <- max(10, ceiling(matrix_profile_size / 100))
  n_work <- floor(matrix_profile_size / per_work)
  idx_work <- list()

  for (i in 1:n_work) {
    idx_st <- (i - 1) * per_work + 1
    if (i == n_work) {
      idx_ed <- matrix_profile_size
    } else {
      idx_ed <- i * per_work
    }
    idx_work[[i]] <- idx_st:idx_ed
  }

  tictac <- Sys.time()

  if (verbose > 0) {
    pb <- utils::txtProgressBar(min = 0, max = n_work, style = 3, width = 80)
  }

  i <- NULL # CRAN NOTE fix
  `%dopar%` <- foreach::`%dopar%` # CRAN NOTE fix

  ## compute the matrix profile
  batch <- foreach(
    i = 1:n_work,
    .verbose = FALSE,
    .inorder = FALSE,
    .multicombine = TRUE,
    .options.snow = opts,
    # .combine = combiner,
    # .errorhandling = 'remove',
    .export = "mass"
  ) %dopar% {
    pro_muls <- matrix(0, length(idx_work[[i]]), n_dim)
    pro_idxs <- matrix(0, length(idx_work[[i]]), n_dim)
    pro_muls_right <- matrix(Inf, length(idx_work[[i]]), n_dim)
    pro_idxs_right <- matrix(-1, length(idx_work[[i]]), n_dim)
    pro_muls_left <- matrix(Inf, length(idx_work[[i]]), n_dim)
    pro_idxs_left <- matrix(-1, length(idx_work[[i]]), n_dim)
    dist_pro <- matrix(0, matrix_profile_size, n_dim)
    last_product <- matrix(0, matrix_profile_size, n_dim)
    drop_value <- matrix(0, 1, n_dim)

    for (j in 1:length(idx_work[[i]])) {
      idx <- idx_work[[i]][j]

      query <- as.matrix(data[idx:(idx + window_size - 1), ])

      if (j == 1) {
        for (k in 1:n_dim) {
          mstomp <- mass(data_fft[, k], query[, k], data_size, window_size, data_mean[, k], data_sd[, k], data_mean[idx, k], data_sd[idx, k])
          dist_pro[, k] <- mstomp$distance_profile
          last_product[, k] <- mstomp$last_product
        }
      } else {
        rep_drop_value <- kronecker(matrix(1, matrix_profile_size - 1, 1), t(drop_value))
        rep_query <- kronecker(matrix(1, matrix_profile_size - 1, 1), t(query[window_size, ]))

        last_product[2:(data_size - window_size + 1), ] <- last_product[1:(data_size - window_size), ] -
          data[1:(data_size - window_size), ] * rep_drop_value +
          data[(window_size + 1):data_size, ] * rep_query


        last_product[1, ] <- first_product[idx, ]

        dist_pro <- 2 * (window_size - (last_product - window_size * data_mean * kronecker(matrix(1, matrix_profile_size, 1), t(data_mean[idx, ]))) /
          (data_sd * kronecker(matrix(1, matrix_profile_size, 1), t(data_sd[idx, ]))))
      }

      dist_pro <- Re(dist_pro)
      drop_value <- query[1, ]

      # apply exclusion zone
      exc_zone_st <- max(1, idx - exclusion_zone)
      exc_zone_ed <- min(matrix_profile_size, idx + exclusion_zone)
      dist_pro[exc_zone_st:exc_zone_ed, ] <- Inf
      dist_pro[data_sd < vars()$eps] <- Inf
      if (skip_location[idx] || any(data_sd[idx, !mask_exc] < vars()$eps)) {
        dist_pro[] <- Inf
      }
      dist_pro[skip_location, ] <- Inf

      # apply dimension "must have" and "exclusion"
      dist_pro[, exc_dim] <- Inf

      if (n_must > 0) {
        mask_must <- rep(FALSE, n_dim)
        mask_must[must_dim] <- TRUE
        dist_pro_must <- dist_pro[, mask_must]
        dist_pro[, mask_must] <- -Inf
      }

      # figure out and store the nearest neighbor
      if (n_dim > 1) {
        dist_pro_sort <- t(apply(dist_pro, 1, sort))
      } # sort by row, left to right
      else {
        dist_pro_sort <- dist_pro
      }

      if (n_must > 0) {
        dist_pro_sort[, 1:n_must] <- dist_pro_must
      }
      dist_pro_cum <- rep(0, matrix_profile_size)
      dist_pro_merg <- rep(0, matrix_profile_size)

      for (k in (max(1, n_must):(n_dim - n_exc))) {
        dist_pro_cum <- dist_pro_cum + dist_pro_sort[, k]
        dist_pro_merg[] <- dist_pro_cum / k
        # left matrix_profile
        if (idx > (exclusion_zone + 1)) {
          min_idx <- which.min(dist_pro_merg[1:(idx - exclusion_zone)])
          min_val <- dist_pro_merg[min_idx]
          pro_muls_left[j, k] <- min_val
          pro_idxs_left[j, k] <- min_idx
        }

        # right matrix_profile
        if (idx < (matrix_profile_size - exclusion_zone)) {
          min_idx <- which.min(dist_pro_merg[(idx + exclusion_zone):matrix_profile_size]) + idx + exclusion_zone - 1
          min_val <- dist_pro_merg[min_idx]
          pro_muls_right[j, k] <- min_val
          pro_idxs_right[j, k] <- min_idx
        }

        # normal matrix_profile
        min_idx <- which.min(dist_pro_merg)
        min_val <- dist_pro_merg[min_idx]
        pro_muls[j, k] <- min_val
        pro_idxs[j, k] <- min_idx
      }
    }

    pro_muls <- sqrt(pro_muls)
    pro_muls_left <- sqrt(pro_muls_left)
    pro_muls_right <- sqrt(pro_muls_right)

    res <- list(pro_muls = pro_muls, pro_idxs = pro_idxs, pro_muls_left = pro_muls_left, pro_idxs_left = pro_idxs_left, pro_muls_right = pro_muls_right, pro_idxs_right = pro_idxs_right, idx = i)

    res
  }

  matrix_profile <- matrix(0, matrix_profile_size, n_dim)
  profile_index <- matrix(0, matrix_profile_size, n_dim)
  left_matrix_profile <- matrix(Inf, matrix_profile_size, n_dim)
  left_profile_index <- matrix(-1, matrix_profile_size, n_dim)
  right_matrix_profile <- matrix(Inf, matrix_profile_size, n_dim)
  right_profile_index <- matrix(-1, matrix_profile_size, n_dim)

  for (i in 1:length(batch)) {
    left_profile_index[idx_work[[batch[[i]]$idx]], ] <- batch[[i]]$pro_idxs_left
    left_matrix_profile[idx_work[[batch[[i]]$idx]], ] <- batch[[i]]$pro_muls_left
    right_profile_index[idx_work[[batch[[i]]$idx]], ] <- batch[[i]]$pro_idxs_right
    right_matrix_profile[idx_work[[batch[[i]]$idx]], ] <- batch[[i]]$pro_muls_right
    profile_index[idx_work[[batch[[i]]$idx]], ] <- batch[[i]]$pro_idxs
    matrix_profile[idx_work[[batch[[i]]$idx]], ] <- batch[[i]]$pro_muls
  }

  ## remove bad k setting in the returned matrix
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
    message(sprintf("\nFinished in %.2f %s", tictac, units(tictac)))
  }

  return(list(
    rmp = right_matrix_profile, rpi = right_profile_index,
    lmp = left_matrix_profile, lpi = left_profile_index,
    mp = matrix_profile, pi = profile_index
  ))
}
