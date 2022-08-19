# [ ]: rlang::check_installed
# [ ]: https://cli.r-lib.org/reference/cli-config.html
# [ ]: https://rlang.r-lib.org/reference/topic-error-call.html


convert_data <- function(data) {
  if (!is.null(data)) {

    try_fetch({
      data_class <- class(data)

      if(is.list(data)) { # tbl_df, data.frame, data.table, list, irts, etc
        if (inherits(data, "irts")) {
          # tseries package
          # data[["time"]]
          # data[["value"]]
        }

      } else if (is.array(data)) { # matrix, array, mts; timeSeries
        if(inherits(data, "timeSeries")) {
          thedata <- timeSeries::getDataPart(data)
          thetime <- timeSeries::getTime(data)
          colnames <- getUnits(data) # same as names(data)
        }

        tspx <- tsp(data)

      } else if (is.ts(data)) {
        tspx <- tsp(data)
      } else if (inherits(data, "xts")) { # matrix, array if unclassed
      # check if irregular time series
        tspx <- xts::xtsAttributes(data)
      } else if (inherits(data, "zoo")) { # matrix, array if unclassed
        tspx <- zoo::index(data)
      }

      # timeSeries, stats::ts, irts (tseries), fts, matrix, data.frame, and zoo::zoo. its

      # if(inherits(data, c("tbl_df", "tbl", "data.frame", "list"))) {
      # } else if(inherits(data, c("matrix", "array"))) {
      # }
    },
     error = function(cnd) {cli::cli_abort("Something went wrong:", parent = cnd)}
     )


    switch()

    if (data_class[[1]] == "numeric") {
      return(invisible(data)) # fast track
    }

    result <- FALSE
    for (class in data_class) {
      if (class %in% c("data.frame", "list")) {
        result <- TRUE
      } else if (class == "array") {
        if (length(dim(data)) > 1) {
          result <- TRUE
        }
      }
    }
    if (result) {
      fn <- rlang::caller_fn() # retrieve the caller function
      if (is.null(fn)) {
        # here this function was called directly
        arg <- rlang::caller_arg(data)
        parenv <- rlang::current_env()
        arg_name <- rlang::fn_fmls_names()[1]
      } else {
        args <- rlang::fn_fmls_names(fn)
        arg_name <- rlang::caller_arg(data)
        i <- which(arg_name == args) + 1
        arg <- rlang::caller_call()[[i]]
        arg <- rlang::as_label(arg)
        parenv <- rlang::caller_env()
      }
      cli::cli_abort("{.arg {arg_name}} = {.var {arg}} cannot be {.cls {class(data)}} and must have a single dimension.", call = parenv)
    }
  }

  return(invisible(as.numeric(data)))
}


#' Clip a value between min and max, fast R implementation
#'
#' For faster yet, use Rcpp::clamp
#'
#' @param x a `vector` or a column `matrix` of `numeric`.
#' @param a a `numeric`, the min value.
#' @param b a `numeric`, the max value.
#'
#' @return Returns the clipped values
#'
#' @keywords internal
#' @noRd
clip <- function(x, a, b) {
  a + (x - a > 0) * (x - a) - (x - b > 0) * (x - b)
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
    cli::cli_abort("`x` and `y` must be numeric vectors or matrices.")
  }
  if (is.vector(x)) {
    dim(x) <- c(1, length(x))
  }
  if (is.vector(y)) {
    dim(y) <- c(1, length(y))
  }
  if (ncol(x) != ncol(y)) {
    cli::cli_abort("`x` and `y` must have the same number of columns.")
  }
  m <- nrow(x)
  n <- nrow(y)
  xy <- x %*% t(y)
  xx <- matrix(rep(apply(x * x, 1, sum), n), m, n, byrow = FALSE)
  yy <- matrix(rep(apply(y * y, 1, sum), m), m, n, byrow = TRUE)
  sqrt(pmax(xx + yy - 2 * xy, 0))
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

# TODO: Revise this algorithm at source

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
#' #
#' # TODO: Refactor
#' #
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
#
# TODO: Refactor
#
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
#
# TODO: Refactor
#
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
#
# TODO: Refactor
#
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
#
# TODO: Refactor
#
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
#
# TODO: Refactor
#
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
#
# TODO: Refactor
#
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
#
# TODO: Refactor
#
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
#
# TODO: Refactor
#
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
#
# TODO: Refactor
#
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
#
# TODO: Refactor
#
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
#
# TODO: Revisit
#
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
#' #
#' # TODO: Refactor
#' #
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
        cli::cli_warn("data size is {nrow(data[[i]])}, but should be {data_size} for this matrix profile.")
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
#' #
#' # TODO: Refactor
#' #
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
#
# TODO: Refactor
#
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
#' #
#' # TODO: Refactor
#' #
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
#' #
#' # TODO: Refactor
#' #
as.matrixprofile <- function(.mp) { # nolint
  if (!("MatrixProfile" %in% class(.mp))) {
    cli::cli_abort("This object cannot be a `MatrixProfile`.")
  }

  class(.mp) <- update_class(class(.mp), "MatrixProfile")

  return(.mp)
}

#' @describeIn as.matrixprofile Cast an object changed by another function back to `MultiMatrixProfile`.
#' @export
#'

as.multimatrixprofile <- function(.mp) { # nolint
  if (!("MultiMatrixProfile" %in% class(.mp))) {
    cli::cli_abort("This object cannot be a `MultiMatrixProfile`.")
  }

  class(.mp) <- update_class(class(.mp), "MultiMatrixProfile")

  return(.mp)
}

#' @describeIn as.matrixprofile Cast an object changed by another function back to `PMP`.
#' @export
#'

as.pmp <- function(.mp) { # nolint
  if (!("PMP" %in% class(.mp))) {
    cli::cli_abort("This object cannot be a `PMP`.")
  }

  class(.mp) <- update_class(class(.mp), "PMP")

  return(.mp)
}

#' @describeIn as.matrixprofile Cast an object changed by another function back to `MultiMatrixProfile`.
#' @export
#'

as.valmod <- function(.mp) { # nolint
  if (!("Valmod" %in% class(.mp))) {
    cli::cli_abort("This object cannot be a `Valmod`.")
  }

  class(.mp) <- update_class(class(.mp), "Valmod")

  return(.mp)
}


#' @describeIn as.matrixprofile Cast an object changed by another function back to `Fluss`.
#' @export

as.fluss <- function(.mp) { # nolint
  if (!("Fluss" %in% class(.mp))) {
    cli::cli_abort("This object cannot be a `Fluss`.")
  }

  class(.mp) <- update_class(class(.mp), "Fluss")

  return(.mp)
}


#' @describeIn as.matrixprofile Cast an object changed by another function back to `Chain`.
#' @export

as.chain <- function(.mp) { # nolint
  if (!("Chain" %in% class(.mp))) {
    cli::cli_abort("This object cannot be a `Chain`.")
  }

  class(.mp) <- update_class(class(.mp), "Chain")

  return(.mp)
}

#' @describeIn as.matrixprofile Cast an object changed by another function back to `Discord`.
#' @export

as.discord <- function(.mp) { # nolint
  if (!("Discord" %in% class(.mp))) {
    cli::cli_abort("This object cannot be a `Discord`.")
  }

  class(.mp) <- update_class(class(.mp), "Discord")

  return(.mp)
}

#' @describeIn as.matrixprofile Cast an object changed by another function back to `Motif`.
#' @export

as.motif <- function(.mp) { # nolint
  if (!("Motif" %in% class(.mp))) {
    cli::cli_abort("This object cannot be a `Motif`.")
  }

  class(.mp) <- update_class(class(.mp), "Motif")

  return(.mp)
}


#' @describeIn as.matrixprofile Cast an object changed by another function back to `MultiMotif`.
#' @export

as.multimotif <- function(.mp) { # nolint
  if (!("MultiMotif" %in% class(.mp))) {
    cli::cli_abort("This object cannot be a `MultiMotif`.")
  }

  class(.mp) <- update_class(class(.mp), "MultiMotif")

  return(.mp)
}


#' @describeIn as.matrixprofile Cast an object changed by another function back to `ArcCount`.
#' @export

as.arccount <- function(.mp) { # nolint
  if (!("ArcCount" %in% class(.mp))) {
    cli::cli_abort("This object cannot be a `ArcCount`.")
  }

  class(.mp) <- update_class(class(.mp), "ArcCount")

  return(.mp)
}

#' @describeIn as.matrixprofile Cast an object changed by another function back to `Salient`.
#' @export

as.salient <- function(.mp) { # nolint
  if (!("Salient" %in% class(.mp))) {
    cli::cli_abort("This object cannot be a `Salient`.")
  }

  class(.mp) <- update_class(class(.mp), "Salient")

  return(.mp)
}
