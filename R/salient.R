#' Retrieve salient subsequences from a dataset
#'
#' In order to allow a meaningful visualization in Multi-Dimensional Space (MDS), this function
#' retrieves the most relevant subsequences using Minimal Description Length (MDL) framework.
#'
#' @details
#' The main purpose of this algorithm is to find subsequences in one time series, but this
#' implementation also covers the experimental effectiveness evaluation with "whole sequence"
#' setting. This means you can input a `matrix` where each column is a sequence and this algorithm
#' will retrieve the most relevant sequences. For this setting, the `exclusion_zone` is ignored, and
#' you need to pre-compute the ordinary euclidean distance matrix. See examples.
#'
#' The `exclusion_zone` is used to avoid trivial matches.
#'
#' `verbose` changes how much information is printed by this function; `0` means nothing,
#' `1` means text, `2` means text and sound.
#'
#' @param data a `vector`, column `matrix` or `data.frame`. If more than one column is provided, see
#'   details.
#' @param matrix_profile a result from STAMP or STOMP algorithms.
#' @param profile_index a result from STAMP or STOMP algorithms.
#' @param window_size an `int` with the size of the sliding window.
#' @param n_bits an `int`. Number of bits for MDL discretization. (Default is `8`).
#' @param n_cand an `int`. number of candidate when picking the subsequence in each iteration.
#'   (Default is `10`).
#' @param exclusion_zone a `numeric`. Size of the exclusion zone, based on `window_size`. (Default
#'   is `1/2`). See details.
#' @param verbose an `int`. See details. (Default is `2`).
#'
#' @return Returns a `list` with `indexes`, a `vector` with the starting position of each
#'   subsequence, `idx_bit_size`, a `vector` with the associated bitsize for each iteration and
#'   `bits` the value used as input on `n_bits`.
#' @references * Yeh CCM, Van Herle H, Keogh E. Matrix profile III: The matrix profile allows
#'   visualization of salient subsequences in massive time series. Proc - IEEE Int Conf Data Mining,
#'   ICDM. 2017;579–88.
#' @references * Hu B, Rakthanmanon T, Hao Y, Evans S, Lonardi S, Keogh E. Discovering the Intrinsic
#'   Cardinality and Dimensionality of Time Series Using MDL. In: 2011 IEEE 11th International
#'   Conference on Data Mining. IEEE; 2011. p. 1086–91.
#' @references Website: <https://sites.google.com/site/salientsubs/>
#' @export
#'
#' @examples
#' \dontrun{
#'   # subsequences setting (main purpose)
#'   salient_subsequences(data, data$mp, data$pi, 30, n_bits = 8, n_cand = 10)
#'
#'   # sequences setting
#'   dist_matrix <- as.matrix(stats::dist(carfull$data))
#'   mp <- matrix(0, nrow(carfull$data), 1)
#'   pi <- matrix(0, nrow(carfull$data), 1)
#'
#'   for (i in 1:nrow(carfull$data)) {
#'     dist_matrix[i, i] <- Inf;
#'     pi[i] <- which.min(dist_matrix[i, ])
#'     mp[i] <- dist_matrix[i, pi[i]]
#'   }
#'
#'   n_bits <- 8
#'
#'   subs <- salient_subsequences(t(carfull$data), mp, pi, 577, n_bits = n_bits, n_cand = 10)
#'   cutoff <- which(diff(subs$idx_bit_size) > 0)[1] - 1
#'   if(cutoff > 0) {
#'     carfull$lab[subs$indexes[1:(cutoff - 1)]]
#'   } else {
#'     message("nothing to do")
#'   }
#' }
#'
salient_subsequences <- function(.mp, data, n_bits = 8, n_cand = 10, exclusion_zone = NULL, verbose = 2) {
  if (!any(class(.mp) %in% c("MatrixProfile"))) {
    stop("Error: First argument must be an object of class `MatrixProfile`.")
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
    stop("Error: Unknown type of data. Must be: matrix, data.frame, vector or list.", call. = FALSE)
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

  matrix_profile <- .mp$mp # keep mp intact

  if (n_dim > 1) {
    max_index_num <- n_dim
  }
  else {
    # get data size
    matrix_profile_size <- nrow(matrix_profile)
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
  uncompressed_bit <- n_bits * .mp$w
  mismatch_bit <- n_bits + log2(.mp$w)

  if (n_dim > 1) {
    idx_bit_size[1] <- uncompressed_bit * n_dim
  } else {
    idx_bit_size[1] <- uncompressed_bit * matrix_profile_size
  }

  if (verbose > 0) {
    pb <- utils::txtProgressBar(min = 0, max = max_index_num, style = 3, width = 80)
    on.exit(close(pb))
  }

  if (verbose > 1) {
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
      hypothesis[hypothesis_count, ] <- discrete_norm(hypothesis[hypothesis_count, ], n_bits, data_max, data_min)
    }

    # get the newest compressiable
    if (compressible_count_old != compressible_count) {
      if (n_dim > 1) {
        compressible[compressible_count, ] <- data[, compressible_idx[compressible_count]]
      } else {
        compressible[compressible_count, ] <- data[compressible_idx[compressible_count]:(compressible_idx[compressible_count] + .mp$w - 1), ]
      }
      compressible[compressible_count, ] <- discrete_norm(compressible[compressible_count, ], n_bits, data_max, data_min)
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

      can <- discrete_norm(can, n_bits, data_max, data_min)

      # test the candiate as hypothesis
      candidate_motif_idx <- .mp$pi[candidate_idx[i]]
      if (n_dim > 1) {
        candidate_motif <- data[, candidate_motif_idx]
      } else {
        candidate_motif <- as.matrix(data[candidate_motif_idx:(candidate_motif_idx + .mp$w - 1), ])
      }
      candidate_motif <- discrete_norm(candidate_motif, n_bits, data_max, data_min)
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

    if (verbose > 0) {
      utils::setTxtProgressBar(pb, indexes_count)
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
  idx_bit_size <- idx_bit_size[1:indexes_count]
  indexes <- indexes[1:indexes_count]

  tictac <- Sys.time() - tictac

  if (verbose > 0) {
    utils::setTxtProgressBar(pb, max_index_num)
    message(sprintf("\nFinished in %.2f %s", tictac, units(tictac)))
  }

  .mp$salient <- list(indexes = indexes, idx_bit_size = idx_bit_size, bits = n_bits)
  class(.mp) <- update_class(class(.mp), "Salient")
  return(.mp)
}

#' Retrieve the index of a number of candidates from the lowest points of a MP
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

# To be ------------------------------------------------------------------------------------------

#' Future function to see MDS
#'
#' @param data original data
#' @param sub_picking picked subsequences
#' @param window_size window size
#'
#' @keywords internal
#' @noRd
#'
#' @return Returns X,Y values for plotting

salient_mds <- function(data, sub_picking, window_size) {
  subs <- list()

  for (i in seq_len(length(sub_picking$indexes))) {
    subs[[i]] <- data[sub_picking$indexes[i]:(sub_picking$indexes[i] + window_size - 1), ]
    subs[[i]] <- (subs[[i]] - mean(subs[[i]])) / std(subs[[i]]) # normalize
  }

  subs <- t(sapply(subs, rbind, simplify = TRUE))
  cmd <- stats::cmdscale(stats::dist(subs), k = 2)

  return(cmd)
}

#' Future function to check performance
#'
#' @param gtruth Ground truth annotation.
#' @param subs Output from `salient_results`.
#' @param window Sliding window size.
#'
#' @return Returns X,Y values for plotting
#'
#' @examples
#' \dontrun{
#'   salient_score(carfull$lab, subs)
#'   salient_score(carsub$labIdx, subssub, carsub$subLen)
#' }
#'
#' @keywords internal
#' @noRd

salient_score <- function(gtruth, subs, window = 0) {
  window <- as.numeric(window)
  best_f <- 0
  best_p <- 0
  best_r <- 0
  best_bit <- 0
  cor_th <- 0.2

  for (i in subs$bits) {
    message("bits: ", i)

    hit_miss <- rep(FALSE, length(subs$indexes))

    for (k in seq_len(length(subs$indexes))) {
      if ((window == 0 && gtruth[subs$indexes[k]] > 0) ||
        (min(abs(subs$indexes[k] - gtruth)) < cor_th * window) # sub
      ) {
        hit_miss[k] <- TRUE
      }
    }

    cutoff <- which(diff(subs$idx_bit_size) > 0)[1] - 1

    if (!is.na(cutoff) && cutoff > 0) {
      hit_miss <- hit_miss[1:cutoff]

      precision <- sum(hit_miss) / length(hit_miss)
      if (window == 0) {
        recall <- sum(hit_miss) / sum(gtruth > 0)
      } else {
        recall <- sum(hit_miss) / length(gtruth)
      }

      f_score <- 2 * precision * recall / (precision + recall)

      if (f_score > best_f) {
        best_p <- precision
        best_r <- recall
        best_f <- f_score
        best_bit <- i
      }

      message("Precision: ", round(precision, 4))
      message("Recall: ", round(recall, 4))
      message("f_score: ", round(f_score, 4))
    } else {
      message("nothing to do")
    }
  }

  message("Best Score: ", round(best_f, 4), " Bits: ", best_bit)

  return(list(precision = best_p, recall = best_r, fscore = f_score, best_bit = best_bit))
}
