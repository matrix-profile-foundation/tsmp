#' Title
#'
#' @param data
#' @param matrix.profile
#' @param profile.index
#' @param window.size
#' @param bits
#' @param n.cand
#' @param exclusion.zone
#'
#' @return
#' @references 1. Yeh CCM, Van Herle H, Keogh E. Matrix profile III: The matrix profile allows visualization of salient subsequences in massive time series. Proc - IEEE Int Conf Data Mining, ICDM. 2017;579–88.
#' @references 2. Hu B, Rakthanmanon T, Hao Y, Evans S, Lonardi S, Keogh E. Discovering the Intrinsic Cardinality and Dimensionality of Time Series Using MDL. In: 2011 IEEE 11th International Conference on Data Mining. IEEE; 2011. p. 1086–91.
#' @references Website: <https://sites.google.com/site/salientsubs/>
#' @export
#'
#' @examples
#' distMat <- as.matrix(dist(carfull$data))
#' matrixProfile <- matrix(0, nrow(carfull$data), 1)
#' profileIndex <- matrix(0, nrow(carfull$data), 1)
#'
#' for (i in 1:nrow(carfull$data)) {
#'   distMat[i, i] <- Inf;
#'   profileIndex[i] <- which.min(distMat[i,])
#'   matrixProfile[i] <- distMat[i,profileIndex[i]]
#' }
bitsave.sub.picking <- function(data, matrix.profile, profile.index, window.size, bits, n.cand, exclusion.zone = 1 / 2, verbose = 2) {

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
    exclusion.zone <- round(window.size * exclusion.zone)
  }

  if (n.dim > 1) {
    max.index.num <- data.size
  }
  else {
    ## get data size
    matrix.profile.size <- nrow(matrix.profile)
    max.index.num <- ceiling(data.size / window.size) ## or round?
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
  uncompressed.bit <- bits * window.size
  mismatch.bit <- bits + log2(window.size)

  if (n.dim > 1) {
    idx.bit.size[1] <- uncompressed.bit * data.size
  } else {
    idx.bit.size[1] <- uncompressed.bit * matrix.profile.size
  }

  tictac <- Sys.time()

  ## iteartivly expend list of hypothesis and compressiable
  while (TRUE) {
    # get the newest hypothesis
    if (hypothesis.count.old != hypothesis.count) {
      if (n.dim > 1) {
        hypothesis[hypothesis.count, ] <- data[hypothesis.idx[hypothesis.count], ]
      } else {
        hypothesis[hypothesis.count, ] <- data[hypothesis.idx[hypothesis.count]:(hypothesis.idx[hypothesis.count] + window.size - 1), ]
      }
      hypothesis[hypothesis.count, ] <- discrete.norm(hypothesis[hypothesis.count, ], bits, data.max, data.min)
    }

    # get the newest compressiable
    if (compressible.count.old != compressible.count) {
      if (n.dim > 1) {
        compressible[compressible.count, ] <- data[compressible.idx[compressible.count], ]
      } else {
        compressible[compressible.count, ] <- data[compressible.idx[compressible.count]:(compressible.idx[compressible.count] + window.size - 1), ]
      }
      compressible[compressible.count, ] <- discrete.norm(compressible[compressible.count, ], bits, data.max, data.min)
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
      for (j in 1:hypothesis.count) {
        if (j == 0) {
          break
        }
        new.descr.length.temp <- get.bitsize(compressible[compressible.count, ] - hypothesis[j, ], mismatch.bit)
        if (new.descr.length.temp < new.descr.length) {
          new.descr.length <- new.descr.length.temp
        }
      }
      compress.cost <- compress.cost + new.descr.length
      hypothesis.cost <- uncompressed.bit * hypothesis.count + compressible.count * log2(hypothesis.count)
      if (n.dim > 1) {
        other.cost <- uncompressed.bit * (data.size - hypothesis.count - compressible.count)
      } else {
        other.cost <- uncompressed.bit * (matrix.profile.size - hypothesis.count - compressible.count)
      }
      idx.bit.size[indexes.count] <- compress.cost + hypothesis.cost + other.cost
    } else if (indexes.count > 1) {
      idx.bit.size[indexes.count] <- idx.bit.size[indexes.count - 1]
    }

    if (indexes.count > 0) {
      message(sprintf("%f", round(idx.bit.size[indexes.count])))
    }

    compressible.count.old <- compressible.count
    hypothesis.count.old <- hypothesis.count

    # stop criteria
    if (n.dim > 1) {
      if (indexes.count >= data.size) {
        break
      }
    }
    else {
      if (indexes.count >= (data.size / window.size)) {
        break
      }

      if (indexes.count > 1 && !(is.infinite(idx.bit.size[indexes.count]) && is.infinite(idx.bit.size[indexes.count - 1]))) {
        if ((idx.bit.size[indexes.count] - idx.bit.size[indexes.count - 1]) > 0) {
          indexes.count <- indexes.count - 1
          break
        }
      }
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
        can <- data[candidate.idx[i], ]
      } else {
        can <- as.matrix(data[candidate.idx[i]:(candidate.idx[i] + window.size - 1), ])
      }

      can <- discrete.norm(can, bits, data.max, data.min)

      # test the candiate as hypothesis
      candidate.motif.idx <- profile.index[candidate.idx[i]]
      if (n.dim > 1) {
        candidate.motif <- data[candidate.motif.idx, ]
      } else {
        candidate.motif <- as.matrix(data[candidate.motif.idx:(candidate.motif.idx + window.size - 1), ])
      }
      candidate.motif <- discrete.norm(candidate.motif, bits, data.max, data.min)
      bitsave.hypothesis <- uncompressed.bit - get.bitsize(candidate.motif - can, mismatch.bit)

      # test the candiate as compressiable
      new.descr.length <- Inf
      for (j in 1:hypothesis.count) {
        if (j == 0) {
          break
        }
        new.descr.length.temp <- get.bitsize(can - hypothesis[j, ], mismatch.bit)
        if (new.descr.length.temp < new.descr.length) {
          new.descr.length <- new.descr.length.temp
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
    indexes[indexes.count] <- candidate.idx[best.candidate]

    if (candidate.bitsave[best.candidate, 2] == 1) {
      hypothesis.count <- hypothesis.count + 1
      hypothesis.idx[hypothesis.count] <- candidate.idx[best.candidate]
      message(sprintf("%f Hypo ", indexes.count))
    } else if (candidate.bitsave[best.candidate, 2] == 2) {
      compressible.count <- compressible.count + 1
      compressible.idx[compressible.count] <- candidate.idx[best.candidate]
      message(sprintf("%f Com ", indexes.count))
    }
  }
  idx.bit.size <- idx.bit.size[1:indexes.count]
  indexes <- indexes[1:indexes.count]

  tictac <- Sys.time() - tictac

  if (verbose > 0) {
    message(sprintf("\nFinished in %.2f %s", tictac, units(tictac)))
  }

  return(list(indexes = indexes, idx.bit.size = idx.bit.size))
}

#' Title
#'
#' @param score
#' @param k
#' @param exclusion.zone
#'
#' @return
#' @export
#'
#' @examples
get.sorted.idx <- function(score, k, exclusion.zone = 0) {
  idx <- sort(score, index.return = TRUE)$ix

  if (exclusion.zone > 0) {
    for (i in 1:length(idx)) {
      if (i > min(k, length(idx))) {
        break
      }
      idx.temp <- idx[(i + 1):length(idx)]
      idx.temp <- idx.temp[abs(idx.temp - idx[i]) >= exclusion.zone]
      idx <- c(idx[1:i], idx.temp)
    }
  }

  idx <- idx[!is.infinite(score[idx])]

  if (k > length(idx)) {
    k <- length(idx)
  }

  idx <- idx[1:k]

  return(idx)
}

#' Title
#'
#' @param x
#' @param mismatch.bit
#'
#' @return
#' @export
#'
#' @examples
get.bitsize <- function(x, mismatch.bit) {
  bit.size <- sum(x != 0) * mismatch.bit

  return(bit.size)
}

#' Title
#'
#' @param data
#' @param window.size
#'
#' @return
#' @export
#'
#' @examples
discrete.norm.pre <- function(data, window.size = 1) {
  if (is.vector(data)) {
    data <- as.matrix(data)
  } else {
    if (ncol(data) > 1) {
      window.size <- 1
    }
  }

  max <- -Inf
  min <- Inf
  for (i in 1:(nrow(data) - window.size + 1)) {
    window <- data[i:(i + window.size - 1), ]
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


### per column?
#' Title
#'
#' @param data
#' @param bits
#' @param max
#' @param min
#'
#' @return
#' @export
#'
#' @examples
discrete.norm <- function(data, bits, max, min) {
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
  data <- round(data * (2^bits - 1)) + 1

  return(data)
}

#' Title
#'
#' @param data
#' @param sub.picking
#' @param window.size
#'
#' @return
#' @export
#'
#' @examples
get.mds <- function(data, sub.picking, window.size) {
  subs <- list()

  for (i in 1:length(sub.picking$indexes)) {
    subs[[i]] <- data[sub.picking$indexes[i]:(sub.picking$indexes[i] + window.size - 1), ]
    subs[[i]] <- (subs[[i]] - mean(subs[[i]])) / std(subs[[i]]) # normalize
  }

  subs <- t(sapply(subs, rbind, simplify = TRUE))
  cmd <- stats::cmdscale(dist(subs), k = 2)

  return(cmd)
}
