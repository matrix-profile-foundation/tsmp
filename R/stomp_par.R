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
#' means text and sound. `exclusion_zone` is used to avoid trivial matches; if a query data is
#' provided (join similarity), this parameter is ignored.
#'
#' @param ... a `matrix` or a `vector`. If a second time series is supplied it will be a join matrix
#'   profile.
#' @param window_size an `int`. Size of the sliding window.
#' @param exclusion_zone a `numeric`. Size of the exclusion zone, based on window size (default is
#'   `1/2`). See details.
#' @param verbose an `int`. See details. (Default is `2`).
#' @param n_workers an `int`. Number of workers for parallel. (Default is `2`).
#'
#' @return Returns the matrix profile `mp` and profile index `pi`. It also returns the left and
#'   right matrix profile `lmp`, `rmp` and profile index `lpi`, `rpi` that may be used to detect
#'   Time Series Chains (Yan Zhu 2018).
#' @export
#'
#' @family matrix profile computations
#' @seealso [stamp()], [stamp_par()]; [mstomp()], [mstomp_par()] for multivariate analysis.
#' @references * Zhu Y, Zimmerman Z, Senobari NS, Yeh CM, Funning G. Matrix Profile II : Exploiting
#'   a Novel Algorithm and GPUs to Break the One Hundred Million Barrier for Time Series Motifs and
#'   Joins. Icdm. 2016 Jan 22;54(1):739–48.
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#' mp <- stomp_par(mp_toy_data$data[1:200,1], window_size = 30, verbose = 0)
#' \dontrun{
#' ref_data <- mp_toy_data$data[,1]
#' query_data <- mp_toy_data$data[,2]
#' # self similarity
#' mp <- stomp_par(ref_data, window_size = 30)
#' # join similarity
#' mp <- stomp_par(ref_data, query_data, window_size = 30)
#' }
stomp_par <- function(..., window_size, exclusion_zone = 1 / 2, verbose = 2, n_workers = 2) {
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

  if (query_size > data_size) {
    stop("Error: Query must be smaller or the same size as reference data.", call. = FALSE)
  }
  if (window_size > query_size / 2) {
    stop("Error: Time series is too short relative to desired window size.", call. = FALSE)
  }
  if (window_size < 4) {
    stop("Error: `window_size` must be at least 4.", call. = FALSE)
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

  ## seperate index into different job
  per_work <- max(10, ceiling(num_queries / 100))
  n_work <- floor(num_queries / per_work)
  idx_work <- list()

  for (i in 1:n_work) {
    idx_st <- (i - 1) * per_work + 1
    if (i == n_work) {
      idx_ed <- num_queries
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
    work_len <- length(idx_work[[i]])
    pro_muls <- matrix(Inf, matrix_profile_size, 1)
    pro_idxs <- matrix(-1, matrix_profile_size, 1)
    pro_muls_right <- matrix(Inf, matrix_profile_size, 1)
    pro_idxs_right <- matrix(-1, matrix_profile_size, 1)
    pro_muls_left <- matrix(Inf, matrix_profile_size, 1)
    pro_idxs_left <- matrix(-1, matrix_profile_size, 1)
    dist_pro <- matrix(0, matrix_profile_size, 1)
    last_product <- matrix(0, matrix_profile_size, 1)
    drop_value <- matrix(0, 1, 1)

    for (j in 1:work_len) {
      # compute the distance profile
      idx_st <- idx_work[[i]][1]
      idx_ed <- idx_work[[i]][work_len]
      idx <- idx_work[[i]][j]

      query_window <- as.matrix(query[idx:(idx + window_size - 1), 1])

      if (j == 1) {
        nn <- mass(
          data_fft, query_window, data_size, window_size, data_mean, data_sd,
          query_mean[idx], query_sd[idx]
        )
        dist_pro[, 1] <- nn$distance_profile
        last_product[, 1] <- nn$last_product
      } else {
        last_product[2:(data_size - window_size + 1), 1] <- last_product[1:(data_size - window_size), 1] -
          data[1:(data_size - window_size), 1] * drop_value +
          data[(window_size + 1):data_size, 1] * query_window[window_size, 1]
        last_product[1, 1] <- first_product[idx, 1]
        dist_pro <- 2 * (window_size - (last_product - window_size * data_mean * query_mean[idx]) /
          (data_sd * query_sd[idx]))
      }

      dist_pro <- Re(sqrt(dist_pro))
      drop_value <- query_window[1, 1]

      # apply exclusion zone
      if (exclusion_zone > 0) {
        exc_st <- max(1, idx - exclusion_zone)
        exc_ed <- min(matrix_profile_size, idx + exclusion_zone)
        dist_pro[exc_st:exc_ed, 1] <- Inf
        dist_pro[data_sd < vars()$eps] <- Inf
        if (skip_location[idx] || any(query_sd[idx] < vars()$eps)) {
          dist_pro[] <- Inf
        }
      }

      # left matrix_profile
      ind <- (dist_pro[idx:matrix_profile_size] < pro_muls_left[idx:matrix_profile_size])
      ind <- c(rep(FALSE, (idx - 1)), ind) # pad left
      pro_muls_left[ind] <- dist_pro[ind]
      pro_idxs_left[which(ind)] <- idx

      # right matrix_profile
      ind <- (dist_pro[1:idx] < pro_muls_right[1:idx])
      ind <- c(ind, rep(FALSE, matrix_profile_size - idx)) # pad right
      pro_muls_right[ind] <- dist_pro[ind]
      pro_idxs_right[which(ind)] <- idx

      # normal matrix_profile
      ind <- (dist_pro < pro_muls)
      pro_muls[ind] <- dist_pro[ind]
      pro_idxs[which(ind)] <- idx
    }

    res <- list(
      pro_muls = pro_muls, pro_idxs = pro_idxs,
      pro_muls_left = pro_muls_left, pro_idxs_left = pro_idxs_left,
      pro_muls_right = pro_muls_right, pro_idxs_right = pro_idxs_right
    )

    res
  }

  matrix_profile <- matrix(Inf, matrix_profile_size, 1)
  profile_index <- matrix(-1, matrix_profile_size, 1)
  left_matrix_profile <- matrix(Inf, matrix_profile_size, 1)
  left_profile_index <- matrix(-1, matrix_profile_size, 1)
  right_matrix_profile <- matrix(Inf, matrix_profile_size, 1)
  right_profile_index <- matrix(-1, matrix_profile_size, 1)

  for (i in 1:length(batch)) {
    ind <- (batch[[i]]$pro_muls < matrix_profile)
    matrix_profile[ind] <- batch[[i]]$pro_muls[ind]
    profile_index[ind] <- batch[[i]]$pro_idxs[ind]

    ind <- (batch[[i]]$pro_muls_left < left_matrix_profile)
    left_matrix_profile[ind] <- batch[[i]]$pro_muls_left[ind]
    left_profile_index[ind] <- batch[[i]]$pro_idxs_left[ind]

    ind <- (batch[[i]]$pro_muls_right < right_matrix_profile)
    right_matrix_profile[ind] <- batch[[i]]$pro_muls_right[ind]
    right_profile_index[ind] <- batch[[i]]$pro_idxs_right[ind]
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
