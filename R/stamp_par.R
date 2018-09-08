#' Anytime univariate STAMP algorithm Parallel version
#'
#' @param n_workers an `int`. Number of workers for parallel. (Default is `2`).
#'
#' @export
#'
#' @describeIn stamp Parallel version.

stamp_par <- function(..., window_size, exclusion_zone = 1 / 2, verbose = 2, s_size = Inf, n_workers = 2) {
  args <- list(...)
  data <- args[[1]]
  if (length(args) > 1) {
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
  matrix_profile_size <- data_size - window_size + 1
  num_queries <- query_size - window_size + 1

  if (window_size > query_size / 2) {
    stop("Error: Time series is too short relative to desired window size.")
  }
  if (window_size < 4) {
    stop("Error: `window_size` must be at least 4.")
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

  matrix_profile <- matrix(Inf, matrix_profile_size, 1)
  profile_index <- matrix(-1, matrix_profile_size, 1)

  if (length(args) > 1) { # no RMP and LMP for joins
    left_matrix_profile <- right_matrix_profile <- NULL
    left_profile_index <- right_profile_index <- NULL
  } else {
    left_matrix_profile <- right_matrix_profile <- matrix_profile
    left_profile_index <- right_profile_index <- profile_index
  }

  ssize <- min(s_size, num_queries)
  order <- sample(1:num_queries, size = ssize)

  tictac <- Sys.time()

  cores <- min(max(2, n_workers), parallel::detectCores())

  if (verbose > 0) {
    message("Warming up parallel with ", cores, " cores.")
  }

  cols <- min(num_queries, 100)

  lines <- 0:(ceiling(ssize / cols) - 1)
  if (verbose > 0) {
    pb <- utils::txtProgressBar(min = 0, max = max(lines), style = 3, width = 80)
  }
  cl <- parallel::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  on.exit(parallel::stopCluster(cl))
  if (verbose > 0) {
    on.exit(close(pb), TRUE)
  }
  if (verbose > 1) {
    on.exit(beep(sounds[[1]]), TRUE)
  }
  # anytime must return the result always
  on.exit(return({
    obj <- list(
      mp = matrix_profile, pi = profile_index,
      rmp = right_matrix_profile, rpi = right_profile_index,
      lmp = left_matrix_profile, lpi = left_profile_index,
      w = window_size,
      ez = ez
    )
    class(obj) <- "MatrixProfile"
    obj
  }), TRUE)

  pre <- mass_pre(data, data_size, query, query_size, window_size = window_size)

  j <- NULL # CRAN NOTE fix
  `%dopar%` <- foreach::`%dopar%` # CRAN NOTE fix

  for (k in lines) {
    batch <- foreach::foreach(
      j = 1:cols,
      # .verbose = FALSE,
      .inorder = FALSE,
      .multicombine = TRUE,
      # .options.snow = opts,
      # .combine = combiner,
      # .errorhandling = 'remove',
      .export = "mass"
    ) %dopar% {
      res <- NULL

      index <- k * cols + j
      if (index <= ssize) {
        i <- order[index]
        nn <- mass(pre$data_fft, query[i:(i + window_size - 1)], data_size, window_size, pre$data_mean, pre$data_sd, pre$query_mean[i], pre$query_sd[i])
        distance_profile <- Re(sqrt(nn$distance_profile))

        # apply exclusion zone
        if (exclusion_zone > 0) {
          exc_st <- max(1, i - exclusion_zone)
          exc_ed <- min(matrix_profile_size, i + exclusion_zone)
          distance_profile[exc_st:exc_ed] <- Inf
          distance_profile[pre$data_sd < vars()$eps] <- Inf
          if (skip_location[i] || any(pre$query_sd[i] < vars()$eps)) {
            distance_profile[] <- Inf
          }
        }

        res <- list(dp = distance_profile, i = i)
      }

      res
    }

    for (i in seq_len(length(batch))) {
      curr <- batch[[i]]$i

      if (!is.null(curr)) {
        if (length(args) == 1) { # no RMP and LMP for joins
          # left matrix_profile
          ind <- (batch[[i]]$dp[curr:matrix_profile_size] < left_matrix_profile[curr:matrix_profile_size])
          ind <- c(rep(FALSE, (curr - 1)), ind) # pad left
          left_matrix_profile[ind] <- batch[[i]]$dp[ind]
          left_profile_index[which(ind)] <- curr

          # right matrix_profile
          ind <- (batch[[i]]$dp[1:curr] < right_matrix_profile[1:curr])
          ind <- c(ind, rep(FALSE, matrix_profile_size - curr)) # pad right
          right_matrix_profile[ind] <- batch[[i]]$dp[ind]
          right_profile_index[which(ind)] <- curr
        }

        # normal matrix_profile
        ind <- (batch[[i]]$dp < matrix_profile)
        matrix_profile[ind] <- batch[[i]]$dp[ind]
        profile_index[which(ind)] <- curr
      }
    }

    if (verbose > 0) {
      utils::setTxtProgressBar(pb, k)
    }
  }

  tictac <- Sys.time() - tictac

  if (verbose > 0) {
    message(sprintf("\nFinished in %.2f %s", tictac, units(tictac)))
  }

  # return() is at on.exit() function
}
