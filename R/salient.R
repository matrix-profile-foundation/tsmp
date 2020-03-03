#' Framework for retrieve salient subsequences from a dataset
#'
#' In order to allow a meaningful visualization in Multi-Dimensional Space (MDS), this function
#' retrieves the most relevant subsequences using Minimal Description Length (MDL) framework.
#'
#' @details
#' `verbose` changes how much information is printed by this function; `0` means nothing,
#' `1` means text, `2` adds the progress bar, `3` adds the finish sound.
#'
#' @param .mp a TSMP object of class `MatrixProfile`.
#' @param data the data used to build the Matrix Profile, if not embedded.
#' @param n_bits an `int` or `vector` of `int`. Number of bits for MDL discretization. (Default is `8`).
#' @param n_cand an `int`. number of candidate when picking the subsequence in each iteration.
#'   (Default is `10`).
#' @param exclusion_zone if a `number` will be used instead of embedded value. (Default is `NULL`).
#' @param verbose an `int`. See details. (Default is `2`).
#'
#' @return Returns the input `.mp` object with a new name `salient`. It contains: `indexes`, a `vector`
#' with the starting position of each subsequence, `idx_bit_size`, a `vector` with the associated
#' bitsize for each iteration and `bits` the value used as input on `n_bits`.
#'
#' @references * Yeh CCM, Van Herle H, Keogh E. Matrix profile III: The matrix profile allows
#'   visualization of salient subsequences in massive time series. Proc - IEEE Int Conf Data Mining,
#'   ICDM. 2017;579-88.
#' @references * Hu B, Rakthanmanon T, Hao Y, Evans S, Lonardi S, Keogh E. Discovering the Intrinsic
#'   Cardinality and Dimensionality of Time Series Using MDL. In: 2011 IEEE 11th International
#'   Conference on Data Mining. IEEE; 2011. p. 1086-91.
#' @references Website: <https://sites.google.com/site/salientsubs/>
#' @export
#'
#' @examples
#' # toy example
#' data <- mp_toy_data$data[, 1]
#' mp <- tsmp(data, window_size = 30, verbose = 0)
#' mps <- salient_subsequences(mp, data, verbose = 0)
#' \dontrun{
#' # full example
#' data <- mp_meat_data$sub$data
#' w <- mp_meat_data$sub$sub_len
#' mp <- tsmp(data, window_size = w, verbose = 2, n_workers = 6)
#' mps <- salient_subsequences(mp, data, n_bits = c(4, 6, 8), verbose = 2)
#' }
#'
salient_subsequences <- function(.mp, data, n_bits = 8, n_cand = 10, exclusion_zone = NULL, verbose = getOption("tsmp.verbose", 2)) {
  if (!("MatrixProfile" %in% class(.mp))) {
    stop("First argument must be an object of class `MatrixProfile`.")
  }

  if ("Valmod" %in% class(.mp)) {
    stop("Function not implemented for objects of class `Valmod`.")
  }

  if (missing(data) && !is.null(.mp$data)) {
    data <- .mp$data[[1]]
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

  if (n_dim > 1) {
    exclusion_zone <- 0
  } else {
    # set trivial match exclusion zone
    if (is.null(exclusion_zone)) {
      exclusion_zone <- round(.mp$w * .mp$ez + vars()$eps)
    } else {
      exclusion_zone <- round(.mp$w * exclusion_zone + vars()$eps)
    }
  }

  if (n_dim > 1) {
    max_index_num <- n_dim
  }
  else {
    # get data size
    matrix_profile_size <- nrow(.mp$mp)
    max_index_num <- round(data_size / .mp$w + vars()$eps)
  }

  # preprocess for discretization
  if (n_dim > 1) {
    minmax <- discrete_norm_pre(data)
  } else {
    minmax <- discrete_norm_pre(data, .mp$w)
  }

  data_min <- minmax$min
  data_max <- minmax$max

  all_idx <- NULL
  all_bit_size <- NULL

  for (b in seq_len(length(n_bits))) {
    matrix_profile <- .mp$mp # keep mp intact

    # initialization various vectors
    indexes <- rep(0, max_index_num)
    idx_bit_size <- rep(0, max_index_num)
    hypothesis_idx <- rep(0, max_index_num)
    compressible_idx <- rep(0, max_index_num)
    hypothesis <- matrix(0, max_index_num, .mp$w)
    compressible <- matrix(0, max_index_num, .mp$w)

    # initialization count related variables
    hypothesis_count <- 0
    compressible_count <- 0
    indexes_count <- 0
    compressible_count_old <- 0
    hypothesis_count_old <- 0

    # initialization bit size related variables
    compress_cost <- 0
    uncompressed_bit <- n_bits[b] * .mp$w
    mismatch_bit <- n_bits[b] + log2(.mp$w)

    if (n_dim > 1) {
      idx_bit_size[1] <- uncompressed_bit * n_dim
    } else {
      idx_bit_size[1] <- uncompressed_bit * matrix_profile_size
    }

    if (verbose > 1) {
      pb <- progress::progress_bar$new(
        format = "Salient [:bar] :percent at :tick_rate it/s, elapsed: :elapsed, eta: :eta",
        clear = FALSE, total = max_index_num, width = 80
      )
    }

    if (verbose > 2) {
      on.exit(beep(sounds[[1]]), TRUE)
    }

    tictac <- Sys.time()

    # iteartivly expend list of hypothesis and compressiable
    while (TRUE) {
      # get the newest hypothesis
      if (hypothesis_count_old != hypothesis_count) {
        if (n_dim > 1) {
          hypothesis[hypothesis_count, ] <- data[, hypothesis_idx[hypothesis_count]]
        } else {
          hypothesis[hypothesis_count, ] <- data[hypothesis_idx[hypothesis_count]:(hypothesis_idx[hypothesis_count] + .mp$w - 1), ]
        }
        hypothesis[hypothesis_count, ] <- discrete_norm(hypothesis[hypothesis_count, ], n_bits[b], data_max, data_min)
      }

      # get the newest compressiable
      if (compressible_count_old != compressible_count) {
        if (n_dim > 1) {
          compressible[compressible_count, ] <- data[, compressible_idx[compressible_count]]
        } else {
          compressible[compressible_count, ] <- data[compressible_idx[compressible_count]:(compressible_idx[compressible_count] + .mp$w - 1), ]
        }
        compressible[compressible_count, ] <- discrete_norm(compressible[compressible_count, ], n_bits[b], data_max, data_min)
      }

      # remove newest hypothesis from the matrix_profile
      if (hypothesis_count_old != hypothesis_count) {
        if (n_dim > 1) {
          matrix_profile[hypothesis_idx[hypothesis_count]] <- Inf
        } else {
          exc_idx_st <- max(1, hypothesis_idx[hypothesis_count] - exclusion_zone)
          exc_idx_ed <- min(matrix_profile_size, hypothesis_idx[hypothesis_count] + exclusion_zone)
          matrix_profile[exc_idx_st:exc_idx_ed] <- Inf
        }
      }

      # remove newest compressiable from the matrix_profile
      if (compressible_count_old != compressible_count) {
        if (n_dim > 1) {
          matrix_profile[compressible_idx[compressible_count]] <- Inf
        } else {
          exc_idx_st <- max(1, compressible_idx[compressible_count] - exclusion_zone)
          exc_idx_ed <- min(matrix_profile_size, compressible_idx[compressible_count] + exclusion_zone)
          matrix_profile[exc_idx_st:exc_idx_ed] <- Inf
        }
      }

      # get current bitsave
      if (compressible_count_old != compressible_count) {
        new_descr_length <- Inf
        if (hypothesis_count > 0) {
          for (j in 1:hypothesis_count) {
            new_descr_length_temp <- get_bitsize(compressible[compressible_count, ] - hypothesis[j, ], mismatch_bit)
            if (new_descr_length_temp < new_descr_length) {
              new_descr_length <- new_descr_length_temp
            }
          }
        }
        compress_cost <- compress_cost + new_descr_length
        hypothesis_cost <- uncompressed_bit * hypothesis_count + compressible_count * log2(hypothesis_count)
        if (n_dim > 1) {
          other_cost <- uncompressed_bit * (n_dim - hypothesis_count - compressible_count)
        } else {
          other_cost <- uncompressed_bit * (matrix_profile_size - hypothesis_count - compressible_count)
        }
        idx_bit_size[indexes_count] <- compress_cost + hypothesis_cost + other_cost
      } else if (indexes_count > 1) {
        idx_bit_size[indexes_count] <- idx_bit_size[indexes_count - 1]
      }

      compressible_count_old <- compressible_count
      hypothesis_count_old <- hypothesis_count

      # stop criteria
      if (indexes_count >= max_index_num) {
        break
      }

      # get candidates
      candidate_idx <- get_sorted_idx(matrix_profile, n_cand, exclusion_zone)
      candidate_idx <- candidate_idx[!is.infinite(matrix_profile[candidate_idx])]
      candidate_n_temp <- length(candidate_idx)

      if (length(candidate_idx) == 0 || any(is.na(candidate_idx))) {
        break
      }

      # testing each candidate
      candidate_bitsave <- matrix(-Inf, candidate_n_temp, 2) # 2nd column (1:hypothesis, 2:compressible)

      for (i in 1:candidate_n_temp) {
        if (n_dim > 1) {
          can <- data[, candidate_idx[i]]
        } else {
          can <- as.matrix(data[candidate_idx[i]:(candidate_idx[i] + .mp$w - 1), ])
        }

        can <- discrete_norm(can, n_bits[b], data_max, data_min)

        # test the candiate as hypothesis
        candidate_motif_idx <- .mp$pi[candidate_idx[i]]
        if (n_dim > 1) {
          candidate_motif <- data[, candidate_motif_idx]
        } else {
          candidate_motif <- as.matrix(data[candidate_motif_idx:(candidate_motif_idx + .mp$w - 1), ])
        }
        candidate_motif <- discrete_norm(candidate_motif, n_bits[b], data_max, data_min)
        bitsave_hypothesis <- uncompressed_bit - get_bitsize(candidate_motif - can, mismatch_bit)

        # test the candiate as compressiable
        new_descr_length <- Inf
        if (hypothesis_count > 0) {
          for (j in 1:hypothesis_count) {
            new_descr_length_temp <- get_bitsize(can - hypothesis[j, ], mismatch_bit)
            if (new_descr_length_temp < new_descr_length) {
              new_descr_length <- new_descr_length_temp
            }
          }
        }
        bitsave_compressed <- uncompressed_bit - new_descr_length

        # if the candidate is better as hypothesis or else
        if (bitsave_hypothesis > bitsave_compressed) {
          candidate_bitsave[i, 1] <- bitsave_hypothesis
          candidate_bitsave[i, 2] <- 1
        } else {
          candidate_bitsave[i, 1] <- bitsave_compressed
          candidate_bitsave[i, 2] <- 2
        }
      }

      # if the candidate is better as hypothesis
      best_candidate <- which.max(candidate_bitsave[, 1])

      if (all(is.infinite(candidate_bitsave[, 2]))) {
        break
      }

      indexes_count <- indexes_count + 1

      if (verbose > 1) {
        pb$tick()
      }

      indexes[indexes_count] <- candidate_idx[best_candidate]

      if (candidate_bitsave[best_candidate, 2] == 1) {
        hypothesis_count <- hypothesis_count + 1
        hypothesis_idx[hypothesis_count] <- candidate_idx[best_candidate]
      } else if (candidate_bitsave[best_candidate, 2] == 2) {
        compressible_count <- compressible_count + 1
        compressible_idx[compressible_count] <- candidate_idx[best_candidate]
      }
    }

    all_bit_size <- cbind(all_bit_size, idx_bit_size[1:indexes_count])
    all_idx <- cbind(all_idx, indexes[1:indexes_count])
  }

  tictac <- Sys.time() - tictac

  if (verbose > 0) {
    message(sprintf("Finished in %.2f %s", tictac, units(tictac)))
  }

  .mp$salient <- list(indexes = as.matrix(all_idx), idx_bit_size = as.matrix(all_bit_size), bits = n_bits)
  class(.mp) <- update_class(class(.mp), "Salient")
  return(.mp)
}


#' Convert salient sequences into MDS space
#'
#' @param .mp a Matrix Profile object.
#' @param data the data used to build the Matrix Profile, if not embedded.
#' @param bit_idx an `int`. The index of `n_bits` used for MDL discretization if more than one was
#' used. (Default is `1`).
#'
#' @return Returns X,Y values for plotting
#'
#' @references * Yeh CCM, Van Herle H, Keogh E. Matrix profile III: The matrix profile allows
#'   visualization of salient subsequences in massive time series. Proc - IEEE Int Conf Data Mining,
#'   ICDM. 2017;579-88.
#' @references * Hu B, Rakthanmanon T, Hao Y, Evans S, Lonardi S, Keogh E. Discovering the Intrinsic
#'   Cardinality and Dimensionality of Time Series Using MDL. In: 2011 IEEE 11th International
#'   Conference on Data Mining. IEEE; 2011. p. 1086-91.
#' @references Website: <https://sites.google.com/site/salientsubs/>
#'
#' @export
#'
#' @examples
#' # toy example
#' data <- mp_toy_data$data[, 1]
#' mp <- tsmp(data, window_size = 30, verbose = 0)
#' mps <- salient_subsequences(mp, verbose = 0)
#' mds_data <- salient_mds(mps)
#' plot(mds_data, main = "Multi dimensional scale")
salient_mds <- function(.mp, data, bit_idx = 1) {
  if (!("Salient" %in% class(.mp))) {
    stop("First argument must be an object of class `Salient`.")
  }

  if (missing(data) && !is.null(.mp$data)) {
    data <- .mp$data[[1]]
  }

  subs <- list()

  for (i in seq_len(nrow(.mp$salient$indexes))) {
    subs[[i]] <- data[.mp$salient$indexes[i, bit_idx]:(.mp$salient$indexes[i, bit_idx] + .mp$w - 1), ]
    subs[[i]] <- (subs[[i]] - mean(subs[[i]])) / std(subs[[i]]) # normalize
  }

  subs <- t(sapply(subs, rbind, simplify = TRUE))
  cmd <- stats::cmdscale(stats::dist(subs), k = 2)
  colnames(cmd) <- c("Dimension 1", "Dimension 2")

  return(cmd)
}

#' Computes the F-Score of salient algorithm.
#'
#' This score function is useful for testing several values of `n_bits` for MDL discretization and
#' checking against a known set of indexes. This increase the probability of better results on
#' relevant subsequence extraction.
#'
#' @param .mp a Matrix Profile object.
#' @param gtruth a `vector` of `integers` with the indexes of relevant subsequences.
#' @param verbose an `int`. (Default is `2`).
#'
#' @return Returns a `list` with `f_score`, `precision`, `recall` and `bits` used in the algorithm.
#'
#' @references * Yeh CCM, Van Herle H, Keogh E. Matrix profile III: The matrix profile allows
#'   visualization of salient subsequences in massive time series. Proc - IEEE Int Conf Data Mining,
#'   ICDM. 2017;579-88.
#' @references * Hu B, Rakthanmanon T, Hao Y, Evans S, Lonardi S, Keogh E. Discovering the Intrinsic
#'   Cardinality and Dimensionality of Time Series Using MDL. In: 2011 IEEE 11th International
#'   Conference on Data Mining. IEEE; 2011. p. 1086-91.
#' @references Website: <https://sites.google.com/site/salientsubs/>
#'
#' @export
#'
#' @examples
#' # toy example
#' data <- mp_toy_data$data[, 1]
#' mp <- tsmp(data, window_size = 30, verbose = 0)
#' mps <- salient_subsequences(mp, n_bits = c(4, 6, 8), verbose = 0)
#' label_idx <- seq(2, 500, by = 110) # fake data
#' salient_score(mps, label_idx, verbose = 0)
salient_score <- function(.mp, gtruth, verbose = getOption("tsmp.verbose", 2)) {
  if (!("Salient" %in% class(.mp))) {
    stop("First argument must be an object of class `Salient`.")
  }

  f_score <- 0
  best_f <- 0
  best_p <- 0
  best_r <- 0
  best_bit <- 0
  cor_th <- 0.2

  for (b in seq_len(length(.mp$salient$bits))) {
    if (verbose > 0) {
      message("Bits: ", .mp$salient$bits[b])
    }

    hit_miss <- rep(FALSE, length(.mp$salient$indexes))

    for (k in seq_len(nrow(.mp$salient$indexes))) {
      if ((.mp$w == 0 && gtruth[.mp$salient$indexes[k, b]] > 0) ||
        (min(abs(.mp$salient$indexes[k, b] - gtruth)) < cor_th * .mp$w) # sub
      ) {
        hit_miss[k] <- TRUE
      }
    }

    cutoff <- which(diff(.mp$salient$idx_bit_size[, b]) > 0)[1] - 1

    if (!is.na(cutoff) && cutoff > 0) {
      hit_miss <- hit_miss[1:cutoff]

      precision <- sum(hit_miss) / length(hit_miss)
      if (.mp$w == 0) {
        recall <- sum(hit_miss) / sum(gtruth > 0)
      } else {
        recall <- sum(hit_miss) / length(gtruth)
      }

      f_score <- 2 * precision * recall / (precision + recall)
      if (is.na(f_score)) {
        f_score <- 0
      }

      if (f_score > best_f) {
        best_p <- precision
        best_r <- recall
        best_f <- f_score
        best_bit <- .mp$salient$bits[b]
      }

      if (verbose > 0) {
        message("Precision: ", round(precision, 4))
        message("Recall: ", round(recall, 4))
        message("F_1 Score: ", round(f_score, 4))
        message("-----------------")
      }
    } else {
      if (verbose > 0) {
        message("Nothing to do")
      }
    }
  }

  if (verbose > 0) {
    message("Best F_1 Score: ", round(best_f, 4), " Bits: ", best_bit)
  }

  return(invisible(list(fscore = f_score, precision = best_p, recall = best_r, bits = best_bit)))
}
