#' Retrieve salient subsequences from a dataset
#'
#' In order to allow a meaningful visualization in Multi-Dimensional Space (MDS), this function
#' retrieves the most relevant subsequences using Minimal Description Length (MDL) framework.
#'
#' The main purpose of this algorithm is to find subsequences in one time series, but this
#' implementation also covers the experimental effectiveness evaluation with "whole sequence"
#' setting. This means you can input a `matrix` where each column is a sequence and this algorithm
#' will retrieve the most relevant sequences. For this setting, the `exclusion.zone` is ignored, and
#' you need to pre-compute the ordinary euclidean distance matrix. See examples.
#'
#' The `exclusion.zone` is used to avoid trivial matches.
#'
#' `verbose` changes how much information is printed by this function; `0` means nothing,
#' `1` means text, `2` means text and sound.
#'
#' @param data a `vector`, column `matrix` or `data.frame`. If more than one column is provided, see
#'   details.
#' @param matrix.profile a result from STAMP or STOMP algorithms.
#' @param profile.index a result from STAMP or STOMP algorithms.
#' @param window.size an `int` with the size of the sliding window.
#' @param n.bits an `int`. Number of bits for MDL discretization. (Default is `8`).
#' @param n.cand an `int`. number of candidate when picking the subsequence in each iteration.
#'   (Default is `10`).
#' @param exclusion.zone a `numeric`. Size of the exclusion zone, based on `window.size`. (Default
#'   is `1/2`). See details.
#' @param verbose an `int`. See details. (Default is `2`).
#'
#' @return Returns a `list` with `indexes`, a `vector` with the starting position of each
#'   subsequence, `idx.bit.size`, a `vector` with the associated bitsize for each iteration and
#'   `bits` the value used as input on `n.bits`.
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
#'   salient.subsequences(data, data$mp, data$pi, 30, n.bits = 8, n.cand = 10)
#'
#'   # sequences setting
#'   dist.matrix <- as.matrix(stats::dist(carfull$data))
#'   mp <- matrix(0, nrow(carfull$data), 1)
#'   pi <- matrix(0, nrow(carfull$data), 1)
#'
#'   for (i in 1:nrow(carfull$data)) {
#'     dist.matrix[i, i] <- Inf;
#'     pi[i] <- which.min(dist.matrix[i, ])
#'     mp[i] <- dist.matrix[i, pi[i]]
#'   }
#'
#'   n.bits <- 8
#'
#'   subs <- salient.subsequences(t(carfull$data), mp, pi, 577, n.bits = n.bits, n.cand = 10)
#'   cutoff <- which(diff(subs$idx.bit.size) > 0)[1] - 1
#'   if(cutoff > 0) {
#'     carfull$lab[subs$indexes[1:(cutoff - 1)]]
#'   } else {
#'     message("nothing to do")
#'   }
#' }
#'
salient.subsequences <- function(data, matrix.profile, profile.index, window.size, n.bits = 8, n.cand = 10, exclusion.zone = 1 / 2, verbose = 2) {

  ## transform data list into matrix
  if (is.matrix(data) || is.data.frame(data)) {
    if (is.data.frame(data)) {
      data <- as.matrix(data)
    } # just to be uniform
    if (ncol(data) > nrow(data)) {
      data <- t(data)
    }
    data.size <- nrow(data)
    n.dim <- ncol(data)
  } else if (is.list(data)) {
    data.size <- length(data[[1]])
    n.dim <- length(data)

    for (i in 1:n.dim) {
      len <- length(data[[i]])
      # Fix TS size with NaN
      if (len < data.size) {
        data[[i]] <- c(data[[i]], rep(NA, data.size - len))
      }
    }
    # transform data into matrix (each column is a TS)
    data <- sapply(data, cbind)
  } else if (is.vector(data)) {
    data.size <- length(data)
    n.dim <- 1
    # transform data into 1-col matrix
    data <- as.matrix(data) # just to be uniform
  } else {
    stop("Unknown type of data. Must be: matrix, data.frame, vector or list")
  }

  if (n.dim > 1) {
    exclusion.zone <- 0
  } else {
    ## set trivial match exclusion zone
    exclusion.zone <- round(window.size * exclusion.zone + vars()$eps)
  }

  if (n.dim > 1) {
    max.index.num <- n.dim
  }
  else {
    ## get data size
    matrix.profile.size <- nrow(matrix.profile)
    max.index.num <- round(data.size / window.size + vars()$eps)
  }

  ## preprocess for discretization
  if (n.dim > 1) {
    minmax <- discrete.norm.pre(data)
  } else {
    minmax <- discrete.norm.pre(data, window.size)
  }

  data.min <- minmax$min
  data.max <- minmax$max

  ## initialization various vectors
  indexes <- rep(0, max.index.num)
  idx.bit.size <- rep(0, max.index.num)
  hypothesis.idx <- rep(0, max.index.num)
  compressible.idx <- rep(0, max.index.num)
  hypothesis <- matrix(0, max.index.num, window.size)
  compressible <- matrix(0, max.index.num, window.size)

  ## initialization count related variables
  hypothesis.count <- 0
  compressible.count <- 0
  indexes.count <- 0
  compressible.count.old <- 0
  hypothesis.count.old <- 0

  ## initialization bit size related variables
  compress.cost <- 0
  uncompressed.bit <- n.bits * window.size
  mismatch.bit <- n.bits + log2(window.size)

  if (n.dim > 1) {
    idx.bit.size[1] <- uncompressed.bit * n.dim
  } else {
    idx.bit.size[1] <- uncompressed.bit * matrix.profile.size
  }

  if (verbose > 0) {
    pb <- utils::txtProgressBar(min = 0, max = max.index.num, style = 3, width = 80)
    on.exit(close(pb))
  }

  if (verbose > 1) {
    on.exit(beep(sounds[[1]]), TRUE)
  }

  tictac <- Sys.time()

  ## iteartivly expend list of hypothesis and compressiable
  while (TRUE) {
    # get the newest hypothesis
    if (hypothesis.count.old != hypothesis.count) {
      if (n.dim > 1) {
        hypothesis[hypothesis.count, ] <- data[, hypothesis.idx[hypothesis.count]]
      } else {
        hypothesis[hypothesis.count, ] <- data[hypothesis.idx[hypothesis.count]:(hypothesis.idx[hypothesis.count] + window.size - 1), ]
      }
      hypothesis[hypothesis.count, ] <- discrete.norm(hypothesis[hypothesis.count, ], n.bits, data.max, data.min)
    }

    # get the newest compressiable
    if (compressible.count.old != compressible.count) {
      if (n.dim > 1) {
        compressible[compressible.count, ] <- data[, compressible.idx[compressible.count]]
      } else {
        compressible[compressible.count, ] <- data[compressible.idx[compressible.count]:(compressible.idx[compressible.count] + window.size - 1), ]
      }
      compressible[compressible.count, ] <- discrete.norm(compressible[compressible.count, ], n.bits, data.max, data.min)
    }

    # remove newest hypothesis from the matrix.profile
    if (hypothesis.count.old != hypothesis.count) {
      if (n.dim > 1) {
        matrix.profile[hypothesis.idx[hypothesis.count]] <- Inf
      } else {
        exc.idx.st <- max(1, hypothesis.idx[hypothesis.count] - exclusion.zone)
        exc.idx.ed <- min(matrix.profile.size, hypothesis.idx[hypothesis.count] + exclusion.zone)
        matrix.profile[exc.idx.st:exc.idx.ed] <- Inf
      }
    }

    # remove newest compressiable from the matrix.profile
    if (compressible.count.old != compressible.count) {
      if (n.dim > 1) {
        matrix.profile[compressible.idx[compressible.count]] <- Inf
      } else {
        exc.idx.st <- max(1, compressible.idx[compressible.count] - exclusion.zone)
        exc.idx.ed <- min(matrix.profile.size, compressible.idx[compressible.count] + exclusion.zone)
        matrix.profile[exc.idx.st:exc.idx.ed] <- Inf
      }
    }

    # get current bitsave
    if (compressible.count.old != compressible.count) {
      new.descr.length <- Inf
      if (hypothesis.count > 0) {
        for (j in 1:hypothesis.count) {
          new.descr.length.temp <- get.bitsize(compressible[compressible.count, ] - hypothesis[j, ], mismatch.bit)
          if (new.descr.length.temp < new.descr.length) {
            new.descr.length <- new.descr.length.temp
          }
        }
      }
      compress.cost <- compress.cost + new.descr.length
      hypothesis.cost <- uncompressed.bit * hypothesis.count + compressible.count * log2(hypothesis.count)
      if (n.dim > 1) {
        other.cost <- uncompressed.bit * (n.dim - hypothesis.count - compressible.count)
      } else {
        other.cost <- uncompressed.bit * (matrix.profile.size - hypothesis.count - compressible.count)
      }
      idx.bit.size[indexes.count] <- compress.cost + hypothesis.cost + other.cost
    } else if (indexes.count > 1) {
      idx.bit.size[indexes.count] <- idx.bit.size[indexes.count - 1]
    }

    compressible.count.old <- compressible.count
    hypothesis.count.old <- hypothesis.count

    # stop criteria
    if (indexes.count >= max.index.num) {
      break
    }

    # get candidates
    candidate.idx <- get.sorted.idx(matrix.profile, n.cand, exclusion.zone)
    candidate.idx <- candidate.idx[!is.infinite(matrix.profile[candidate.idx])]
    candidate.n.temp <- length(candidate.idx)

    if (length(candidate.idx) == 0 || any(is.na(candidate.idx))) {
      break
    }

    # testing each candidate
    candidate.bitsave <- matrix(-Inf, candidate.n.temp, 2) # 2nd column (1:hypothesis, 2:compressible)

    for (i in 1:candidate.n.temp) {
      if (n.dim > 1) {
        can <- data[, candidate.idx[i]]
      } else {
        can <- as.matrix(data[candidate.idx[i]:(candidate.idx[i] + window.size - 1), ])
      }

      can <- discrete.norm(can, n.bits, data.max, data.min)

      # test the candiate as hypothesis
      candidate.motif.idx <- profile.index[candidate.idx[i]]
      if (n.dim > 1) {
        candidate.motif <- data[, candidate.motif.idx]
      } else {
        candidate.motif <- as.matrix(data[candidate.motif.idx:(candidate.motif.idx + window.size - 1), ])
      }
      candidate.motif <- discrete.norm(candidate.motif, n.bits, data.max, data.min)
      bitsave.hypothesis <- uncompressed.bit - get.bitsize(candidate.motif - can, mismatch.bit)

      # test the candiate as compressiable
      new.descr.length <- Inf
      if (hypothesis.count > 0) {
        for (j in 1:hypothesis.count) {
          new.descr.length.temp <- get.bitsize(can - hypothesis[j, ], mismatch.bit)
          if (new.descr.length.temp < new.descr.length) {
            new.descr.length <- new.descr.length.temp
          }
        }
      }
      bitsave.compressed <- uncompressed.bit - new.descr.length

      # if the candidate is better as hypothesis or else
      if (bitsave.hypothesis > bitsave.compressed) {
        candidate.bitsave[i, 1] <- bitsave.hypothesis
        candidate.bitsave[i, 2] <- 1
      } else {
        candidate.bitsave[i, 1] <- bitsave.compressed
        candidate.bitsave[i, 2] <- 2
      }
    }

    # if the candidate is better as hypothesis
    best.candidate <- which.max(candidate.bitsave[, 1])

    if (all(is.infinite(candidate.bitsave[, 2]))) {
      break
    }

    indexes.count <- indexes.count + 1

    if (verbose > 0) {
      utils::setTxtProgressBar(pb, indexes.count)
    }

    indexes[indexes.count] <- candidate.idx[best.candidate]

    if (candidate.bitsave[best.candidate, 2] == 1) {
      hypothesis.count <- hypothesis.count + 1
      hypothesis.idx[hypothesis.count] <- candidate.idx[best.candidate]
    } else if (candidate.bitsave[best.candidate, 2] == 2) {
      compressible.count <- compressible.count + 1
      compressible.idx[compressible.count] <- candidate.idx[best.candidate]
    }
  }
  idx.bit.size <- idx.bit.size[1:indexes.count]
  indexes <- indexes[1:indexes.count]

  tictac <- Sys.time() - tictac

  if (verbose > 0) {
    utils::setTxtProgressBar(pb, max.index.num)
    message(sprintf("\nFinished in %.2f %s", tictac, units(tictac)))
  }

  return(list(indexes = indexes, idx.bit.size = idx.bit.size, bits = n.bits))
}

#' Retrieve the index of a number of candidates from the lowest points of a MP
#'
#' @param matrix.profile the matrix profile
#' @param n.cand number of candidates to extract
#' @param exclusion.zone exclusion zone for extracting candidates (in absolute values)
#'
#' @return Returns the indexes of candidates
#'
#' @keywords internal
#'
get.sorted.idx <- function(matrix.profile, n.cand, exclusion.zone = 0) {
  idx <- sort(matrix.profile, index.return = TRUE)$ix

  if (exclusion.zone > 0) {
    for (i in 1:length(idx)) {
      if (i > min(n.cand, length(idx))) {
        break
      }
      idx.temp <- idx[(i + 1):length(idx)]
      idx.temp <- idx.temp[abs(idx.temp - idx[i]) >= exclusion.zone]
      idx <- c(idx[1:i], idx.temp)
    }
  }

  idx <- idx[!is.infinite(matrix.profile[idx])]

  if (n.cand > length(idx)) {
    n.cand <- length(idx)
  }

  idx <- idx[1:n.cand]

  return(idx)
}

#' Reduced description length
#'
#' @param x the difference between two time series (reference and candidate for compression)
#' @param mismatch.bit sum of n.bits and log2(window.size)
#'
#' @return Returns the bit.size cost of compressing the time series
#' @keywords internal

get.bitsize <- function(x, mismatch.bit) {
  bit.size <- sum(x != 0) * mismatch.bit

  return(bit.size)
}

#' Precompute the max and min value for the discrete normalization
#'
#' @param data input time series
#' @param window.size sliding window size
#'
#' @return Returns a list with the max and min value
#' @keywords internal
#'
discrete.norm.pre <- function(data, window.size = 1) {
  if (is.vector(data)) {
    data <- as.matrix(data)
  }

  if (ncol(data) > 1) {
    len <- ncol(data)
  } else {
    len <- nrow(data) - window.size + 1
  }

  max <- -Inf
  min <- Inf
  for (i in 1:len) {
    if (ncol(data) > 1) {
      window <- data[, i]
    } else {
      window <- data[i:(i + window.size - 1), ]
    }
    window.mean <- mean(window)
    window.sd <- std(window)
    if (window.sd == 0) {
      window <- (window - window.mean)
    } else {
      window <- (window - window.mean) / window.sd
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
#' @param n.bits Number of bits for MDL discretization.
#' @param max Precomputed max from `discrete.norm.pre`.
#' @param min Precomputed min from `discrete.norm.pre`.
#'
#' @return Returns the data after discrete normalization.
#' @keywords internal

discrete.norm <- function(data, n.bits, max, min) {
  # normalize magnitude
  data.mean <- mean(data)
  data.sd <- std(data)

  if (data.sd == 0) {
    data <- (data - data.mean)
  } else {
    data <- (data - data.mean) / data.sd
  }

  data <- (data - min) / (max - min)

  # discretization
  data <- round(data * (2^n.bits - 1) + vars()$eps) + 1

  return(data)
}

#' Future function to see MDS
#'
#' @param data original data
#' @param sub.picking picked subsequences
#' @param window.size window size
#'
#' @keywords internal
#'
#' @return Returns X,Y values for plotting

salient.mds <- function(data, sub.picking, window.size) {
  subs <- list()

  for (i in 1:length(sub.picking$indexes)) {
    subs[[i]] <- data[sub.picking$indexes[i]:(sub.picking$indexes[i] + window.size - 1), ]
    subs[[i]] <- (subs[[i]] - mean(subs[[i]])) / std(subs[[i]]) # normalize
  }

  subs <- t(sapply(subs, rbind, simplify = TRUE))
  cmd <- stats::cmdscale(stats::dist(subs), k = 2)

  return(cmd)
}

#' Future function to check performance
#'
#' @param gtruth Ground truth annotation.
#' @param subs Output from `salient.results`.
#' @param window Sliding window size.
#'
#' @return Returns X,Y values for plotting
#'
#' @examples
#' \dontrun{
#'   salient.score(carfull$lab, subs)
#'   salient.score(carsub$labIdx, subssub, carsub$subLen)
#' }
#'
#' @keywords internal

salient.score <- function(gtruth, subs, window = 0) {
  window <- as.numeric(window)
  best.f <- 0
  best.p <- 0
  best.r <- 0
  best.bit <- 0
  cor.th <- 0.2

  for (i in subs$bits) {
    message("bits: ", i)

    hit.miss <- rep(FALSE, length(subs$indexes))

    for (k in 1:length(subs$indexes)) {
      if ((window == 0 && gtruth[subs$indexes[k]] > 0) ||
        (min(abs(subs$indexes[k] - gtruth)) < cor.th * window) # sub
      ) {
        hit.miss[k] <- TRUE
      }
    }

    cutoff <- which(diff(subs$idx.bit.size) > 0)[1] - 1

    if (!is.na(cutoff) && cutoff > 0) {
      hit.miss <- hit.miss[1:cutoff]

      precision <- sum(hit.miss) / length(hit.miss)
      if (window == 0) {
        recall <- sum(hit.miss) / sum(gtruth > 0)
      } else {
        recall <- sum(hit.miss) / length(gtruth)
      }

      f.score <- 2 * precision * recall / (precision + recall)

      if (f.score > best.f) {
        best.p <- precision
        best.r <- recall
        best.f <- f.score
        best.bit <- i
      }

      message("Precision: ", round(precision, 4))
      message("Recall: ", round(recall, 4))
      message("f.score: ", round(f.score, 4))
    } else {
      message("nothing to do")
    }
  }

  message("Best Score: ", round(best.f, 4), " Bits: ", best.bit)

  return(list(precision = best.p, recall = best.r, fscore = f.score, best.bit = best.bit))
}
