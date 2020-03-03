# Window functions ----------------------------------------------------------------------------

# Supports:
#               NA/NaN  -Inf/Inf  Edge  Rcpp
# movmin          Unk     Unk      No   Yes
# movmax          Unk     Unk      No   Yes
# fast_movsd      No      Unk      No   No
# fast_movavg     No      Unk      No   No
# fast_avg_sd     No      Unk      No   No

#' Fast implementation of moving standard deviation
#'
#' This function does not handle NA values
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#' @param window_size moving sd window size
#' @param rcpp a `logical`. Uses rcpp implementation.
#'
#' @return Returns a `vector` with the moving standard deviation
#' @export
#'
#' @examples
#' data_sd <- fast_movsd(mp_toy_data$data[, 1], mp_toy_data$sub_len)
fast_movsd <- function(data, window_size, rcpp = FALSE) {
  if (window_size < 2) {
    stop("'window_size' must be at least 2.")
  }

  # Rcpp is slower
  if (rcpp) {
    return(fast_movsd_rcpp(data, window_size))
  }

  # Improve the numerical analysis by subtracting off the series mean
  # this has no effect on the standard deviation.
  data <- data - mean(data)

  data_sum <- cumsum(c(sum(data[1:window_size]), diff(data, window_size)))
  data_mean <- data_sum / window_size

  data2 <- data^2
  data2_sum <- cumsum(c(sum(data2[1:window_size]), diff(data2, window_size)))
  data_sd2 <- (data2_sum / window_size) - (data_mean^2) # variance
  data_sd <- sqrt(data_sd2)

  return(data_sd)
}

#' Fast implementation of moving average
#'
#' This function does not handle NA values
#'
#' @inheritParams fast_movsd
#'
#' @return Returns a `vector` with the moving average
#' @export
#'
#' @examples
#' data_avg <- fast_movavg(mp_toy_data$data[, 1], mp_toy_data$sub_len)
fast_movavg <- function(data, window_size) {
  if (window_size < 2) {
    stop("'window_size' must be at least 2.")
  }

  return(cumsum(c(sum(data[1:window_size]), diff(data, window_size))) / window_size)
}

#' Converts euclidean distances into correlation values
#'
#' @param x a `vector` or a column `matrix` of `numeric`.
#' @param w the window size
#'
#' @return Returns the converted values
#'
#' @keywords internal
#' @noRd
ed_corr <- function(x, w) {
  (2 * w - x^2) / (2 * w)
}

#' Converts correlation values into euclidean distances
#'
#' @inheritParams ed_corr
#'
#' @return Returns the converted values
#'
#' @keywords internal
#' @noRd
corr_ed <- function(x, w) {
  sqrt(2 * w * (1 - ifelse(x > 1, 1, x)))
}

#' Fast implementation of moving average and moving standard deviation
#'
#' This function does not handle NA values
#'
#' @inheritParams fast_movsd
#'
#' @return Returns a `list` with `avg` and `sd` `vector`s
#' @export

fast_avg_sd <- function(data, window_size, rcpp = FALSE) {
  if (window_size < 2) {
    stop("'window_size' must be at least 2.")
  }

  # Rcpp is slower
  if (rcpp) {
    return(fast_avg_sd_rcpp(data, window_size))
  }

  mov_sum <- cumsum(c(sum(data[1:window_size]), diff(data, window_size)))
  data2 <- data^2
  mov2_sum <- cumsum(c(sum(data2[1:window_size]), diff(data2, window_size)))
  mov_mean <- mov_sum / window_size


  # Improve the numerical analysis by subtracting off the series mean
  # this has no effect on the standard deviation.
  dmean <- mean(data)
  data <- data - dmean

  data_sum <- cumsum(c(sum(data[1:window_size]), diff(data, window_size)))
  data_mean <- data_sum / window_size
  data2 <- data^2
  data2_sum <- cumsum(c(sum(data2[1:window_size]), diff(data2, window_size)))
  data_sd2 <- (data2_sum / window_size) - (data_mean^2) # variance
  data_sd2[data_sd2 < 0] <- 0
  data_sd <- sqrt(data_sd2) # std deviation
  data_sig <- sqrt(1 / (data_sd2 * window_size))

  return(list(avg = mov_mean, sd = data_sd, sig = data_sig, sum = mov_sum, sqrsum = mov2_sum))
}

# DO NOT Handles NA's
fast_muinvn <- function(data, window_size) {
  if (window_size < 2) {
    stop("'window_size' must be at least 2.")
  }

  data_sum <- cumsum(c(sum(data[1:window_size]), diff(data, window_size)))
  data_mean <- data_sum / window_size
  data2 <- data^2
  data2_sum <- cumsum(c(sum(data2[1:window_size]), diff(data2, window_size)))
  data_dp <- 1 / sqrt(data2_sum - data_mean^2 * window_size)

  return(list(avg = data_mean, isd = data_dp))
}

# Handles NA's
old_fast_movsd <- function(data, window_size) {

  # length of the time series
  data_size <- length(data)

  if (window_size < 2) {
    stop("'window_size' must be at least 2.")
  }

  if (data_size < window_size) {
    stop("'window_size' is too large for this series.")
  }

  # Improve the numerical analysis by subtracting off the series mean
  # this has no effect on the standard deviation.
  data <- data - mean(data, na.rm = TRUE)

  # scale the data to have unit variance too. will put that
  # scale factor back into the result at the end
  data_sd <- std(data, na.rm = TRUE)
  data <- data / data_sd

  # we will need the squared elements
  data_sqr <- data^2

  b <- rep(1, window_size)
  s <- sqrt((stats::filter(data_sqr, b, sides = 1) - (stats::filter(data, b, sides = 1)^2) * (1 / window_size)) / (window_size - 1))

  # restore the scale factor that was used before to normalize the data
  s <- s * data_sd
  s <- Re(s)
  s <- s * sqrt((window_size - 1) / window_size)

  return(s[window_size:data_size])
}

# Handles NA's
old_fast_movavg <- function(data, window_size) {
  if (window_size < 2) {
    stop("'window_size' must be at least 2.")
  }

  data_mean <- stats::filter(data, rep(1 / window_size, window_size), sides = 2)
  return(data_mean[!is.na(data_mean)])
}

# DO NOT Handles NA's
# catastrophic cancellation
old_fast_avg_sd <- function(data, window_size) {
  if (window_size < 2) {
    stop("'window_size' must be at least 2.")
  }

  data_len <- length(data)

  data[(data_len + 1):(window_size + data_len)] <- 0

  data_cum <- cumsum(data)
  data2_cum <- cumsum(data^2)
  data2_sum <- data2_cum[window_size:data_len] - c(0, data2_cum[1:(data_len - window_size)])
  data_sum <- data_cum[window_size:data_len] - c(0, data_cum[1:(data_len - window_size)])

  data_mean <- data_sum / window_size

  data_sd2 <- (data2_sum / window_size) - (data_mean^2)
  data_sd2 <- pmax(data_sd2, 0)
  data_sd <- sqrt(data_sd2)

  return(list(avg = data_mean, sd = data_sd, sum = data_sum, sqrsum = data2_sum))
}


# Math functions ------------------------------------------------------------------------------
#
# Supports:
#               NA/NaN  -Inf/Inf
# std             Unk      Unk
# mode            Unk      Unk
# znorm           Unk      Unk
# diff2           Unk      Unk
# bubble_up       Unk      Unk
# paa             Unk      Unk
# ipaa            Unk      Unk
# min_mp_idx      Unk      Unk


#' Population SD, as R always calculate with n-1 (sample), here we fix it
#'
#' @inheritParams fast_movsd
#'
#' @return Returns the corrected standard deviation from sample to population
#' @keywords internal
#' @noRd
#'
std <- function(data, na.rm = FALSE, rcpp = TRUE) {

  # Rcpp is faster
  if (rcpp) {
    return(std_rcpp(data, na.rm))
  }

  sdx <- stats::sd(data, na.rm)

  if (is.na(sdx)) {
    return(NA)
  }

  return(sqrt((length(data) - 1) / length(data)) * sdx)
}

#' Calculates the mode of a vector
#'
#' @param x
#'
#' @return the mode
#' @keywords internal
#' @noRd

mode <- function(x, rcpp = FALSE) {

  # Rcpp is not faster
  if (rcpp) {
    return(mode_rcpp(x))
  }

  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' Normalizes data for mean Zero and Standard Deviation One
#'
#' @inheritParams fast_movsd
#'
#' @return Returns the normalized data
#' @keywords internal
#' @noRd
#'
znorm <- function(data, rcpp = TRUE) {
  # Rcpp is faster
  if (rcpp) {
    return(as.matrix(znorm_rcpp(data)))
  }

  data_mean <- mean(data)
  data_dev <- std(data)

  if (is.na(data_dev) || data_dev <= 0.01) {
    return(data - data_mean)
  }
  else {
    (data - data_mean) / (data_dev)
  }
}

#' Normalizes data to be between min and max
#'
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#' @param min the minimum value
#' @param max the maximum value
#'
#' @return Returns the normalized data
#' @keywords internal
#' @noRd
#'
normalize <- function(data, min = 0, max = 1) {
  min_val <- min(data, na.rm = TRUE)
  max_val <- max(data, na.rm = TRUE)

  a <- (max - min) / (max_val - min_val)
  b <- max - a * max_val
  data <- a * data + b

  data[data < min] <- min
  data[data > max] <- max

  return(data)
}

#' Distance between two matrices
#'
#' Computes the Euclidean distance between rows of two matrices.
#'
#' @param x a `matrix`.
#' @param y a `matrix`.
#'
#' @return Returns a `matrix` of size m x n if x is of size m x k and y is of size n x k.
#' @keywords internal
#' @noRd

diff2 <- function(x, y) {
  # Rcpp ?
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("`x` and `y` must be numeric vectors or matrices.")
  }
  if (is.vector(x)) {
    dim(x) <- c(1, length(x))
  }
  if (is.vector(y)) {
    dim(y) <- c(1, length(y))
  }
  if (ncol(x) != ncol(y)) {
    stop("`x` and `y` must have the same number of columns.")
  }
  m <- nrow(x)
  n <- nrow(y)
  xy <- x %*% t(y)
  xx <- matrix(rep(apply(x * x, 1, sum), n), m, n, byrow = FALSE)
  yy <- matrix(rep(apply(y * y, 1, sum), m), m, n, byrow = TRUE)
  sqrt(pmax(xx + yy - 2 * xy, 0))
}

#' Binary Split algorithm
#'
#' Creates a vector with the indexes of binary split.
#'
#' @param n size of the vector
#'
#' @return Returns a `vector` with the binary split indexes
#' @keywords internal
#' @noRd

binary_split <- function(n, rcpp = TRUE) {
  if (rcpp) {
    return(binary_split_rcpp(as.integer(n)))
  }

  if (n < 2) {
    return(1)
  }

  split <- function(lb, ub, m) {
    if (lb == m) {
      l <- NULL
      r <- c(m + 1, ub)
    } else if (ub == m) {
      l <- c(lb, m - 1)
      r <- NULL
    } else {
      l <- c(lb, m - 1)
      r <- c(m + 1, ub)
    }

    return(list(l = l, r = r))
  }

  idxs <- vector(mode = "numeric", length = n)
  intervals <- list()

  idxs[1] <- 1 # We always begin by explore the first integer
  intervals[[1]] <- c(2, n) # After exploring the first integer, we begin splitting the interval 2:n
  i <- 2

  while (length(intervals) > 0) {
    lb <- intervals[[1]][1]
    ub <- intervals[[1]][2]
    mid <- floor((lb + ub) / 2)
    intervals[[1]] <- NULL

    idxs[i] <- mid
    i <- i + 1

    if (lb == ub) {
      next
    } else {
      lr <- split(lb, ub, mid)
      if (!is.null(lr$l)) {
        intervals[[length(intervals) + 1]] <- lr$l
      }
      if (!is.null(lr$r)) {
        intervals[[length(intervals) + 1]] <- lr$r
      }
    }
  }
  return(idxs)
}

#' Bubble up algorithm
#'
#' Bubble up algorithm.
#'
#' @param data a vector of values
#' @param len size of data
#'
#' @return Doesnt return. Not used for now
#' @keywords internal
#' @noRd

bubble_up <- function(data, len) {
  # Rcpp ?
  pos <- len

  while (pos > 0 && data[pos / 2] < data[pos]) {
    # &&
    #    heap->heap_[pos / 2].distance >= 0 &&
    #   heap->heap_[pos].distance >= 0) {
    t <- data[pos]
    data[pos] <- data[pos / 2]
    data[pos / 2] <- t
    pos <- pos / 2
  }
}

#' Piecewise Aggregate Approximation of time series
#'
#' @param data time series
#' @param p factor of PAA reduction (2 == half of size)
#'
#' @return PAA result
#' @keywords internal
#' @noRd

paa <- function(data, p) {
  # Rcpp ?
  paa_data <- as.vector(data)
  len <- length(paa_data)

  p <- round(abs(p))
  paa_size <- len / p

  if (len == paa_size) {
    return(data)
  } else {
    if (len %% paa_size == 0) {
      res <- colMeans(matrix(paa_data, nrow = len %/% paa_size, byrow = FALSE))
    } else {
      stop("Invalid paa_size")
    }
  }

  if (is.matrix(data)) {
    return(as.matrix(res))
  } else {
    return(res)
  }
}

#' Resample data to the original size, with interpolation
#'
#' @param data time series
#' @param p factor of PAA reduction (2 == half of size)
#'
#' @keywords internal
#' @noRd

ipaa <- function(data, p) {
  # Rcpp ?
  if (is.null(data)) {
    return(NULL)
  }

  paa_data <- as.vector(data)
  paa_size <- length(paa_data)
  size <- paa_size * p

  res <- rep.int(NA, size)

  j <- 1
  for (i in seq_len(size)) {
    if (((i - 1) %% p) == 0) {
      res[i] <- data[j]
      j <- j + 1
    } else {
      res[i] <- res[i - 1]
    }
  }

  if (is.matrix(data)) {
    return(as.matrix(res))
  } else {
    return(res)
  }
}

#' Get index of the minimum value from a matrix profile and its nearest neighbor
#'
#' @param .mp a `MatrixProfile` object.
#' @param n_dim number of dimensions of the matrix profile
#' @param valid check for valid numbers
#'
#' @return returns a `matrix` with two columns: the minimum and the nearest neighbor
#' @export
#'
#' @examples
#' w <- 50
#' data <- mp_gait_data
#' mp <- tsmp(data, window_size = w, exclusion_zone = 1 / 4, verbose = 0)
#' min_val <- min_mp_idx(mp)
min_mp_idx <- function(.mp, n_dim = NULL, valid = TRUE) {
  if (!is.null(n_dim)) {
    .mp$mp <- .mp$mp[, n_dim, drop = FALSE]
    .mp$pi <- .mp$pi[, n_dim, drop = FALSE]
  }

  n_dim <- ncol(.mp$mp)
  mp_size <- nrow(.mp$mp)
  min <- apply(.mp$mp, 2, which.min) # support for multidimensional matrix profile

  if (any(min == 1) && any(is.infinite(.mp$mp[1, (min == 1)]))) {
    return(NA)
  }

  nn_min <- NULL

  for (i in seq_len(n_dim)) {
    nn_min <- c(nn_min, .mp$pi[min[i], i])
  }

  if (valid) {
    if (all(nn_min > 0 & nn_min <= mp_size) &&
      all(!is.infinite(diag(.mp$mp[nn_min, ], names = FALSE)))) {
      return(cbind(min, nn_min, deparse.level = 0))
    }

    for (i in seq_len(n_dim)) {
      .mp$mp[min[i], i] <- Inf
    }

    stop <- FALSE
    while (!stop) {
      min <- apply(.mp$mp, 2, which.min)

      if (any(min == 1) && any(is.infinite(.mp$mp[1, (min == 1)]))) {
        stop <- TRUE
      } else {
        nn_min <- NULL

        for (i in seq_len(n_dim)) {
          nn_min <- c(nn_min, .mp$pi[min[i], i])
        }

        if (all(nn_min > 0 & nn_min <= mp_size) &&
          all(!is.infinite(diag(.mp$mp[nn_min, ], names = FALSE)))) {
          return(cbind(min, nn_min, deparse.level = 0))
        } else {
          for (i in seq_len(n_dim)) {
            .mp$mp[min[i], i] <- Inf
          }
        }
      }
    }

    return(NA)
  } else {
    return(cbind(min, nn_min, deparse.level = 0))
  }
}

# SDTS Aux functions -----------------------------------------------------------------------------


#' Computes the golden section for individual candidates
#'
#' @param dist_pro the candidate distance profile
#' @param label a vector with the data bool annotation
#' @param pos_st a vector with the starting points of label
#' @param pos_ed a vector with the ending points of label
#' @param beta a number that balance the F-Score. Beta > 1 towards recall, < towards precision
#' @param window_size an integer with the sliding window size
#'
#' @return Returns the best threshold and its F-Score
#'
#' @keywords internal
#' @noRd
#'
golden_section <- function(dist_pro, label, pos_st, pos_ed, beta, window_size) {
  # Rcpp ?
  golden_ratio <- (1 + sqrt(5)) / 2
  a_thold <- min(dist_pro)
  b_thold <- max(dist_pro[!is.infinite(dist_pro)])
  c_thold <- b_thold - (b_thold - a_thold) / golden_ratio
  d_thold <- a_thold + (b_thold - a_thold) / golden_ratio
  tol <- max((b_thold - a_thold) * 0.001, 0.0001)

  if (anyNA(c(c_thold, d_thold, tol))) {
    return(list(thold = NA, score = 0))
  }

  while (abs(c_thold - d_thold) > tol) {
    c_score <- compute_f_meas(label, pos_st, pos_ed, dist_pro, c_thold, window_size, beta)
    d_score <- compute_f_meas(label, pos_st, pos_ed, dist_pro, d_thold, window_size, beta)

    if (c_score$f_meas > d_score$f_meas) {
      b_thold <- d_thold
    } else {
      a_thold <- c_thold
    }

    c_thold <- b_thold - (b_thold - a_thold) / golden_ratio
    d_thold <- a_thold + (b_thold - a_thold) / golden_ratio
  }
  thold <- (a_thold + b_thold) * 0.5
  score <- compute_f_meas(label, pos_st, pos_ed, dist_pro, thold, window_size, beta)

  return(list(thold = thold, score = score$f_meas))
}

#' Computes the golden section for combined candidates

#' @param dist_pro the candidate distance profile
#'
#' @param thold a number with the threshold used to calculate the F-Score
#' @param label a vector with the data bool annotation
#' @param pos_st a vector with the starting points of label
#' @param pos_ed a vector with the ending points of label
#' @param beta a number that balance the F-Score. Beta > 1 towards recall, < towards precision
#' @param window_size an integer with the sliding window size
#' @param fit_idx an integer with the index of the current threshold
#'
#' @return Returns the best threshold and its F-Score
#'
#' @keywords internal
#' @noRd

golden_section_2 <- function(dist_pro, thold, label, pos_st, pos_ed, beta, window_size, fit_idx) {
  # Rcpp ?
  golden_ratio <- (1 + sqrt(5)) / 2
  a_thold <- min(dist_pro[[fit_idx]])
  b_thold <- max(dist_pro[[fit_idx]][!is.infinite(dist_pro[[fit_idx]])])
  c_thold <- b_thold - (b_thold - a_thold) / golden_ratio
  d_thold <- a_thold + (b_thold - a_thold) / golden_ratio
  tol <- max((b_thold - a_thold) * 0.001, 0.0001)

  if (anyNA(c(c_thold, d_thold, tol))) {
    return(list(thold = NA, score = 0))
  }

  while (abs(c_thold - d_thold) > tol) {
    c_thold_combined <- thold
    d_thold_combined <- thold
    c_thold_combined[fit_idx] <- c_thold
    d_thold_combined[fit_idx] <- d_thold

    c_score <- compute_f_meas(label, pos_st, pos_ed, dist_pro, c_thold_combined, window_size, beta)
    d_score <- compute_f_meas(label, pos_st, pos_ed, dist_pro, d_thold_combined, window_size, beta)

    if (c_score$f_meas > d_score$f_meas) {
      b_thold <- d_thold
    } else {
      a_thold <- c_thold
    }

    c_thold <- b_thold - (b_thold - a_thold) / golden_ratio
    d_thold <- a_thold + (b_thold - a_thold) / golden_ratio
  }
  thold[fit_idx] <- (a_thold + b_thold) * 0.5
  score <- compute_f_meas(label, pos_st, pos_ed, dist_pro, thold, window_size, beta)

  # beta = 2;   emphacise recall
  # beta = 0.5; emphacise precision
  return(list(thold = thold, score = score$f_meas))
}

#' Computes de F-Score
#'
#' @param label a vector with the data bool annotation
#' @param pos_st a vector with the starting points of label
#' @param pos_ed a vector with the ending points of label
#' @param dist_pro the distance profile of the data
#' @param thold a number with the threshold used to compute the prediction
#' @param window_size an integer with the sliding window size
#' @param beta a number that balance the F-Score. Beta > 1 towards recall, < towards precision
#'
#' @return Returns the F-Score, precision and recall values
#'
#' @keywords internal
#' @noRd

compute_f_meas <- function(label, pos_st, pos_ed, dist_pro, thold, window_size, beta) {
  # Rcpp ?
  # generate annotation curve for each pattern
  if (is.list(dist_pro)) {
    anno_st <- list()
    n_pat <- length(dist_pro)

    for (i in 1:n_pat) {
      annor <- dist_pro[[i]] - thold[i]
      annor[annor > 0] <- 0
      annor[annor < 0] <- -1
      annor <- -annor
      anno_st[[i]] <- which(diff(c(0, annor, 0)) == 1) + 1
      anno_st[[i]] <- anno_st[[i]] - 1
    }

    anno_st <- unlist(anno_st)
    anno_st <- sort(anno_st)

    i <- 1
    while (TRUE) {
      if (i >= length(anno_st)) {
        break
      }

      first_part <- anno_st[1:i]
      second_part <- anno_st[(i + 1):length(anno_st)]
      bad_st <- abs(second_part - anno_st[i]) < window_size

      second_part <- second_part[!bad_st]
      anno_st <- c(first_part, second_part)
      i <- i + 1
    }

    anno_ed <- anno_st + window_size - 1
  } else {
    anno <- dist_pro - thold
    anno[anno > 0] <- 0
    anno[anno < 0] <- -1
    anno <- -anno

    anno_st <- which(diff(c(0, anno, 0)) == 1) + 1
    anno_ed <- anno_st + window_size - 1
    anno_st <- anno_st - 1
    anno_ed <- anno_ed - 1
  }

  anno <- rep(FALSE, length(label))

  for (i in seq_len(length(anno_st))) {
    anno[anno_st[i]:anno_ed[i]] <- 1
  }

  is_tp <- rep(FALSE, length(anno_st))

  for (i in seq_len(length(anno_st))) {
    if (anno_ed[i] > length(label)) {
      anno_ed[i] <- length(label)
    }
    if (sum(label[anno_st[i]:anno_ed[i]]) > (0.8 * window_size)) {
      is_tp[i] <- TRUE
    }
  }
  tp_pre <- sum(is_tp)

  is_tp <- rep(FALSE, length(pos_st))
  for (i in seq_len(length(pos_st))) {
    if (sum(anno[pos_st[i]:pos_ed[i]]) > (0.8 * window_size)) {
      is_tp[i] <- TRUE
    }
  }
  tp_rec <- sum(is_tp)

  pre <- tp_pre / length(anno_st)
  rec <- tp_rec / length(pos_st)

  f_meas <- (1 + beta^2) * (pre * rec) / ((beta^2) * pre + rec)
  if (is.na(f_meas)) {
    f_meas <- 0
  }
  return(list(f_meas = f_meas, pre = pre, rec = rec))
}


# Salient Aux functions --------------------------------------------------------------------------

#' Retrieve the index of a number of candidates from the lowest points of a Matrix Profile
#'
#' @param matrix_profile the matrix profile
#' @param n_cand number of candidates to extract
#' @param exclusion_zone exclusion zone for extracting candidates (in absolute values)
#'
#' @return Returns the indexes of candidates
#'
#' @keywords internal
#' @noRd
#'
get_sorted_idx <- function(matrix_profile, n_cand, exclusion_zone = 0) {
  # Rcpp ?
  idx <- sort(matrix_profile, index.return = TRUE)$ix

  if (exclusion_zone > 0) {
    for (i in seq_len(length(idx))) {
      if (i > min(n_cand, length(idx))) {
        break
      }
      idx_temp <- idx[(i + 1):length(idx)]
      idx_temp <- idx_temp[abs(idx_temp - idx[i]) >= exclusion_zone]
      idx <- c(idx[1:i], idx_temp)
    }
  }

  idx <- idx[!is.infinite(matrix_profile[idx])]

  if (n_cand > length(idx)) {
    n_cand <- length(idx)
  }

  idx <- idx[1:n_cand]

  return(idx)
}

#' Reduced description length
#'
#' @param x the difference between two time series (reference and candidate for compression)
#' @param mismatch_bit sum of n_bits and log2(window_size)
#'
#' @return Returns the bit_size cost of compressing the time series
#' @keywords internal
#' @noRd

get_bitsize <- function(x, mismatch_bit) {
  bit_size <- sum(x != 0) * mismatch_bit

  return(bit_size)
}

#' Precompute the max and min value for the discrete normalization
#'
#' @param data input time series
#' @param window_size sliding window size
#'
#' @return Returns a list with the max and min value
#' @keywords internal
#' @noRd
#'
discrete_norm_pre <- function(data, window_size = 1) {
  if (is.vector(data)) {
    data <- as.matrix(data)
  }

  if (ncol(data) > 1) {
    len <- ncol(data)
  } else {
    len <- nrow(data) - window_size + 1
  }

  max <- -Inf
  min <- Inf
  for (i in 1:len) {
    if (ncol(data) > 1) {
      window <- data[, i]
    } else {
      window <- data[i:(i + window_size - 1), ]
    }
    window_mean <- mean(window)
    window_sd <- std(window)
    if (window_sd == 0) {
      window <- (window - window_mean)
    } else {
      window <- (window - window_mean) / window_sd
    }

    if (max(window) > max) {
      max <- max(window)
    }
    if (min(window) < min) {
      min <- min(window)
    }
  }
  return(list(max = max, min = min))
}


#' Discrete normalization
#'
#' @param data Input time series.
#' @param n_bits Number of bits for MDL discretization.
#' @param max Precomputed max from `discrete_norm_pre`.
#' @param min Precomputed min from `discrete_norm_pre`.
#'
#' @return Returns the data after discrete normalization.
#' @keywords internal
#' @noRd

discrete_norm <- function(data, n_bits, max, min) {
  # normalize magnitude
  data_mean <- mean(data)
  data_sd <- std(data)

  if (data_sd == 0) {
    data <- (data - data_mean)
  } else {
    data <- (data - data_mean) / data_sd
  }

  data <- (data - min) / (max - min)

  # discretization
  data <- round(data * (2^n_bits - 1) + vars()$eps) + 1

  return(data)
}

# AV Aux functions --------------------------------------------------------------------------------

#' Counts number of zero-crossings
#'
#' Count the number of zero-crossings from the supplied time-domain input vector. A simple method is
#' applied here that can be easily ported to a real-time system that would minimize the number of
#' if-else conditionals.
#'
#' @param data is the input time-domain signal (one dimensional).
#'
#' @return Returns the amount of zero-crossings in the input signal.
#' @author sparafucile17 06/27/04
#' @references <https://www.dsprelated.com/showcode/179.php>
#' @keywords internal
#' @noRd
#'
zero_crossings <- function(data) {
  # initial value
  count <- 0

  data <- as.matrix(data)

  # error checks
  if (length(data) == 1) {
    stop("Input signal must have more than one element.")
  }

  if ((ncol(data) != 1) && (nrow(data) != 1)) {
    stop("Input must be one-dimensional.")
  }

  # force signal to be a vector oriented in the same direction
  data <- as.vector(data)

  num_samples <- length(data)

  for (i in 2:num_samples) {
    # Any time you multiply to adjacent values that have a sign difference
    # the result will always be negative.  When the signs are identical,
    # the product will always be positive.
    if ((data[i] * data[i - 1]) < 0) {
      count <- count + 1
    }
  }

  return(count)
}

#' Normalizes data between Zero and One
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#'
#' @return Returns the normalized data.
#' @keywords internal
#' @noRd
#'
zero_one_norm <- function(data) {
  data <- round(data, 10)

  data <- data - min(data[!is.infinite(data) & !is.na(data)])
  data <- data / max(data[!is.infinite(data) & !is.na(data)])

  return(data)
}

#' Computes the complexity index of the data
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#'
#' @return Returns the complexity index of the data provided (normally a subset)
#' @keywords internal
#' @noRd
#'
complexity <- function(data) {
  return(sqrt(sum(diff(data)^2)))
}

# Find Motif Ung Aux Functions -------------------------------------------------------------

#' Computes the bits saved using basic compactation algorithm
#'
#' @param motif_1 reference motif
#' @param motif_2 difference motif
#' @param n_dim data dimensions
#' @param n_bit bits for discretization
#'
#' @return Returns the amount of bits saved
#' @keywords internal
#' @noRd
#'
get_bit_save <- function(motif_1, motif_2, n_dim, n_bit) {
  if (is.vector(motif_1)) {
    motif_1 <- as.matrix(motif_1)
  }

  if (is.vector(motif_2)) {
    motif_2 <- as.matrix(motif_2)
  }

  tot_dim <- dim(motif_1)[2]
  window_size <- dim(motif_1)[1]
  split_pt <- get_desc_split_pt(n_bit)
  disc_1 <- discretization(motif_1, split_pt)
  disc_2 <- discretization(motif_2, split_pt)

  dim_id <- sort(apply(abs(disc_1 - disc_2), 2, sum), index.return = TRUE)$ix
  dim_id <- dim_id[1:n_dim]
  motif_diff <- disc_1[, dim_id] - disc_2[, dim_id]
  n_val <- length(unique(as.vector(motif_diff)))

  bit_sz <- n_bit * (tot_dim * window_size * 2 - n_dim * window_size)
  bit_sz <- bit_sz + n_dim * window_size * log2(n_val) + n_val * n_bit

  return(list(bit_sz = bit_sz, dim_id = dim_id))
}

#' Discretize a time series using split points.
#'
#' @param motif a `matrix` with the motif.
#' @param split_pt split points for discretization.
#'
#' @return Returns the discretized time series.
#' @keywords internal
#' @noRd
#'
discretization <- function(motif, split_pt) {
  if (is.vector(motif)) {
    motif <- as.matrix(motif)
  }

  dimmotif <- dim(motif)

  for (i in 1:dimmotif[2]) {
    motif[, i] <- (motif[, i] - mean(motif[, i])) / std(motif[, i])
  }

  disc <- matrix(0, dimmotif[1], dimmotif[2])

  for (i in seq_len(length(split_pt))) {
    disc[motif < split_pt[i] & disc == 0] <- i
  }

  disc[disc == 0] <- length(split_pt) + 1

  return(disc)
}

#' Get the split points for discretization
#'
#' @param n_bit number of bits for discretization.
#'
#' @return Returns the split points.
#' @keywords internal
#' @noRd
#'
get_desc_split_pt <- function(n_bit) {
  split_pt <- stats::qnorm((1:((2^n_bit) - 1)) / (2^n_bit), 0, 1)
  return(split_pt)
}

# Global constants ---------------------------------------------------------------------------

#' Global constants
#'
#' @return Returns a `list` with the global constants
#' @keywords internal
#' @noRd
#'
vars <- function() {
  eps <- .Machine$double.eps^0.5
  kmode <- 0.6311142 # mode is ((a-1) / (a*b-1))^(1/a) ==> 0.6311142

  return(list(eps = eps, kmode = kmode))
}

# Misc -------------------------------------------------------------------------------------------
#' Set/changes the data included in TSMP object.
#'
#' This may be useful if you want to include the data lately or remove the included data (set as `NULL`).
#'
#' @param .mp a TSMP object.
#' @param data a `matrix` (for one series) or a `list` of matrices (for two series).
#'
#' @return Returns silently the original TSMP object with changed data.
#' @export
#'
#' @examples
#' mp <- tsmp(mp_toy_data$data[1:200, 1], window_size = 30, verbose = 0)
#' mp <- set_data(mp, NULL)
set_data <- function(.mp, data) {
  if (!is.null(data)) {
    if (!is.list(data)) {
      data <- list(data)
    }

    data <- lapply(data, as.matrix)

    # data should be this size
    data_size <- (nrow(.mp$mp) + .mp$w - 1)

    for (i in seq_len(length(data))) {
      if (nrow(data[[i]]) != data_size) {
        warning("Warning: data size is ", nrow(data[[i]]), ", but should be ", data_size, " for this matrix profile.")
      }
    }
  }

  .mp$data <- data

  invisible(.mp)
}

#' Get the data included in a TSMP object, if any.
#'
#' @param .mp a TSMP object.
#'
#' @return Returns the data as `matrix`. If there is more than one series, returns a `list`.
#' @export
#'
#' @examples
#' mp <- tsmp(mp_toy_data$data[1:200, 1], window_size = 30, verbose = 0)
#' get_data(mp)
get_data <- function(.mp) {
  if (length(.mp$data) == 1) {
    return(invisible(as.matrix(.mp$data[[1]])))
  } else {
    return(invisible(as.matrix(.mp$data)))
  }
}

#' Add class on front or move it to front if already exists
#'
#' @param classes result from `class()`
#' @param new_class string with the new class to add or update (put on front).
#'
#' @return Returns a vector with classes names.
#' @keywords internal
#' @noRd
#'

update_class <- function(classes, new_class) {
  classes <- classes[!(classes == new_class)]
  classes <- c(new_class, classes)

  return(classes)
}

#' Remove a `TSMP` class from an object
#'
#' @param x a `TSMP` object
#' @param class `character` string with the class name
#'
#' @return the object without the class
#' @export
#'
#' @examples
#' w <- 50
#' data <- mp_gait_data
#' mp <- tsmp(data, window_size = w, exclusion_zone = 1 / 4, verbose = 0)
#' mp <- find_chains(mp)
#' # Remove the "Chain" class information
#' mp <- remove_class(mp, "Chain")
remove_class <- function(x, class) {
  switch(class,
    "Chain" = {
      x$chain <- NULL
    },
    "Discord" = {
      x$discord <- NULL
    },
    "Motif" = {
      x$motif <- NULL
    },
    "MultiMotif" = {
      x$motif <- NULL
    },
    "AnnotationVector" = {
      x$av <- NULL
    },
    "ArcCount" = {
      x$cac <- NULL
    },
    "Fluss" = {
      x$fluss <- NULL
    },
    "Salient" = {
      x$salient <- NULL
    }
  )

  class(x) <- class(x)[class(x) != class]

  return(x)
}

#' Convert a TSMP object into another if possible
#'
#' The base Classes are `MatrixProfile` and `MultiMatrixProfile`, but as other functions are used,
#' classes are pushed behind, since the last output normally is the most significant. If you want,
#' for example, to plot the Matrix Profile from a `Fluss` object, you may use `as.matrixprofile()`
#' to cast it back.
#'
#' @param .mp a TSMP object.
#'
#' @return Returns the object with the new class, if possible.
#'
#' @describeIn as.matrixprofile Cast an object changed by another function back to `MatrixProfile`.
#' @export
#' @examples
#'
#' w <- 50
#' data <- mp_gait_data
#' mp <- tsmp(data, window_size = w, exclusion_zone = 1 / 4, verbose = 0)
#' mp <- find_motif(mp)
#' class(mp) # first class will be "Motif"
#'
#' plot(mp) # plots a motif plot
#'
#' plot(as.matrixprofile(mp)) # plots a matrix profile plot
as.matrixprofile <- function(.mp) {
  if (!("MatrixProfile" %in% class(.mp))) {
    stop("This object cannot be a `MatrixProfile`.")
  }

  class(.mp) <- update_class(class(.mp), "MatrixProfile")

  return(.mp)
}

#' @describeIn as.matrixprofile Cast an object changed by another function back to `MultiMatrixProfile`.
#' @export
#'

as.multimatrixprofile <- function(.mp) {
  if (!("MultiMatrixProfile" %in% class(.mp))) {
    stop("This object cannot be a `MultiMatrixProfile`.")
  }

  class(.mp) <- update_class(class(.mp), "MultiMatrixProfile")

  return(.mp)
}

#' @describeIn as.matrixprofile Cast an object changed by another function back to `PMP`.
#' @export
#'

as.pmp <- function(.mp) {
  if (!("PMP" %in% class(.mp))) {
    stop("This object cannot be a `PMP`.")
  }

  class(.mp) <- update_class(class(.mp), "PMP")

  return(.mp)
}

#' @describeIn as.matrixprofile Cast an object changed by another function back to `MultiMatrixProfile`.
#' @export
#'

as.valmod <- function(.mp) {
  if (!("Valmod" %in% class(.mp))) {
    stop("This object cannot be a `Valmod`.")
  }

  class(.mp) <- update_class(class(.mp), "Valmod")

  return(.mp)
}


#' @describeIn as.matrixprofile Cast an object changed by another function back to `Fluss`.
#' @export

as.fluss <- function(.mp) {
  if (!("Fluss" %in% class(.mp))) {
    stop("This object cannot be a `Fluss`.")
  }

  class(.mp) <- update_class(class(.mp), "Fluss")

  return(.mp)
}


#' @describeIn as.matrixprofile Cast an object changed by another function back to `Chain`.
#' @export

as.chain <- function(.mp) {
  if (!("Chain" %in% class(.mp))) {
    stop("This object cannot be a `Chain`.")
  }

  class(.mp) <- update_class(class(.mp), "Chain")

  return(.mp)
}

#' @describeIn as.matrixprofile Cast an object changed by another function back to `Discord`.
#' @export

as.discord <- function(.mp) {
  if (!("Discord" %in% class(.mp))) {
    stop("This object cannot be a `Discord`.")
  }

  class(.mp) <- update_class(class(.mp), "Discord")

  return(.mp)
}

#' @describeIn as.matrixprofile Cast an object changed by another function back to `Motif`.
#' @export

as.motif <- function(.mp) {
  if (!("Motif" %in% class(.mp))) {
    stop("This object cannot be a `Motif`.")
  }

  class(.mp) <- update_class(class(.mp), "Motif")

  return(.mp)
}


#' @describeIn as.matrixprofile Cast an object changed by another function back to `MultiMotif`.
#' @export

as.multimotif <- function(.mp) {
  if (!("MultiMotif" %in% class(.mp))) {
    stop("This object cannot be a `MultiMotif`.")
  }

  class(.mp) <- update_class(class(.mp), "MultiMotif")

  return(.mp)
}


#' @describeIn as.matrixprofile Cast an object changed by another function back to `ArcCount`.
#' @export

as.arccount <- function(.mp) {
  if (!("ArcCount" %in% class(.mp))) {
    stop("This object cannot be a `ArcCount`.")
  }

  class(.mp) <- update_class(class(.mp), "ArcCount")

  return(.mp)
}

#' @describeIn as.matrixprofile Cast an object changed by another function back to `Salient`.
#' @export

as.salient <- function(.mp) {
  if (!("Salient" %in% class(.mp))) {
    stop("This object cannot be a `Salient`.")
  }

  class(.mp) <- update_class(class(.mp), "Salient")

  return(.mp)
}

#' Play sound with `audio`
#'
#' @param data sound data provided by this package
#'
#' @keywords internal
#' @noRd
#'
beep <- function(data) {
  if (!(is.null(audio::audio.drivers()) || nrow(audio::audio.drivers()) == 0)) {
    tryCatch(
      {
        audio::play(data)
      },
      error = function(cond) {
        message("Warning: Failed to play audio alert")
        message(cond)
      },
      warning = function(cond) {
        message("Warning: Something went wrong playing audio alert")
        message(cond)
      }
    )
  }
  Sys.sleep(1)
  invisible()
}
