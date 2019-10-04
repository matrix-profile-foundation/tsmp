#' Title
#'
#' @param data
#' @param window_sizes
#'
#'   The work closest in spirit to ours is VALMOD. The idea of VALMOD is to compute the MP for
#'   the shortest length of interest, then use the information gleaned from it to guide a search
#'   through longer subsequence lengths, exploiting lower bounds to prune off some calculations.
#'   This idea works well for the first few of the longer subsequence lengths, but the lower bounds
#'   progressively weaken, making the pruning ineffective. Thus, in the five case studies they
#'   presented, the mean value of U/L was just 1.24. In contrast, consider that our termite example
#'   in Fig. 15 has a U/L ratio of 240, more than two orders of magnitude larger. Thus, VALMOD is
#'   perhaps best seen as finding motifs with some tolerance for a slightly (~25%) too short
#'   user-specified query length, rather than a true "motif-of-all-lengths" algorithm. Also note
#'   that apart from the shortest length, VALMOD only gives some information for the other lengths,
#'   unlike pmp, which contains exact distances for all subsequences of all lengths.
#'
#'
#' @return
#' @export
#'
#' @examples
skimp <- function(data, window_sizes = seq.int(1, seq_end, jump), plot = TRUE) {
  #   skimp generates the Pan Matrix Profile (pmp), a data structure which contains the Matrix
  #   Profile information for all subsequence lengths in a range. If not ran to completion, the
  #   instruction pmp(idx) yields the approximation as a matrix
  #
  #   INPUTS data - A nonempty, row or column vector, consisting of real values. window_sizes   - A
  #   positive and increasing sequence of integers which enumerate the subsequence lengths used to
  #   create the Pan Matrix Profile
  #
  #   OUTPUTS pmp - Pan Matrix Profile approximation. Each row contains a single run of the Matrix
  #   Profile for some subsequence length. If ran to completion this variable stores the exact Pan
  #   Matrix Profile; otherwise, pmp(idx) yields the approximation idx - The permutation of the rows
  #   in the pmp which correspond to the pmp approximation
  #

  # assert(max(INPUTS[['window_sizes']]) < length(data) / 2, '[skimp] Error: The maximum subsequence
  # length explored must be less than half the size of the dataset.')
  #
  # INITIALIZE

  # Determine the order in which we will explore the subsequence lengths used to run the matrix profile
  split_idx <- binary_split(length(window_sizes))

  # Stores the intermediate values of each Matrix Profile where pmp(i,:) = MatrixProfile with
  # subsequence length INPUTS[['window_sizes']](split_idx(i))
  pmp <- matrix(Inf, nrow = length(window_sizes), ncol = length(data))

  # Stores row indices for image (Image is not a copy of the pmp but rather a permutation (with
  # repetition) of its rows
  idx <- rep(-1, length(window_sizes))

  ## BEGIN PROCESSING

  # Runs matrix profile on each subsequence length in INPUTS[['window_sizes']](split_idx), stores the
  # results in pmp, and generates a frame for the animation
  for (i in 1:length(split_idx)) {

    # Get the current subsequence length to run the Matrix Profile on
    window_size <- window_sizes[split_idx[i]]

    # Run Matrix Profile
    mp <- mpx(data, window_size)$mp
    # to do ed
    pmp[split_idx[i], 1:length(mp)] <- mp

    # Propogate changes to idx
    # TO-DO: Use a binary search to find last index, update all values in range. This MAY result in a time save
    j <- split_idx[i]
    while (j <= length(split_idx) && idx[j] != j) {
      idx[j] <- split_idx[i]
      j <- j + 1
    }

    if (plot) {
      test <- pmp[rev(idx), ]
      depth <- 256
      test <- ceiling(test * depth) / depth
      test[test > 1] <- 1
      grid::grid.raster(test, height = 0.3, width = 1, interpolate = FALSE)
    }

    message("i: ", i, "/", length(split_idx))
  }

  pmp <- pmp[rev(idx), ]

  return(list(pmp = pmp, idx = idx))
} # function skimp

#' Title
#'
#' @param n
#'
#' @return
#' @export
#'
#' @examples
binary_split <- function(n, rcpp = TRUE) {
  if (rcpp) {
    return(binary_split_rcpp(as.integer(n)))
  }

  # binary_split(Length)
  # Preallocate memory to store the traversed indices and the remaining intervals
  #
  #
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

#' Title
#'
#' @param lb
#' @param ub
#' @param m
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param a
#' @param minlag
#' @param w
#'
#' @return
#' @export
#'
#' @examples
mpx <- function(a, w, rcpp = TRUE) {
  minlag <- floor(w / 4)

  if (rcpp) {
    return(mpx_rcpp(as.numeric(a), as.integer(w), as.integer(minlag)))
  }

  # matrix profile using cross correlation,
  n <- length(a)

  msd <- fast_avg_sd(a, w)
  mu <- msd$avg
  sig <- msd$sig

  # differentials have 0 as their first entry. This simplifies index
  # calculations slightly and allows us to avoid special "first line"
  # handling.
  df <- c(0, (1 / 2) * (a[(1 + w):n] - a[1:(n - w)]))
  dg <- c(0, (a[(1 + w):n] - mu[2:(n - w + 1)]) + (a[1:(n - w)] - mu[1:(n - w)]))
  diagmax <- length(a) - w + 1
  mp <- rep(-1, n - w + 1)
  mpi <- rep(NA, n - w + 1)

  seq_diag <- (minlag + 1):diagmax
  seq_diag <- sample(seq_diag, size = length(seq_diag))

  # anytime must return the result always
  on.exit(return({
    # to do ed
    mp[mp > 1] <- 1.0
    mp <- sqrt(2 * w * (1 - mp))
    list(mp = mp, mpi = mpi)
  }), TRUE)

  for (diag in seq_diag) {
    c <- sum((a[diag:(diag + w - 1)] - mu[diag]) * (a[1:w] - mu[1]))

    offset <- 1:(n - w - diag + 2)
    off_diag <- (offset + diag - 1)
    d <- df[offset] * dg[off_diag] + df[off_diag] * dg[offset]

    d <- tail(diffinv(d, xi = c), -1)

    d_cmp <- d * sig[offset] * sig[off_diag]

    mask <- d_cmp > mp[offset]
    mp[c(mask, rep(FALSE, diag - 1))] <- d_cmp[mask]
    mpi[c(mask, rep(FALSE, diag - 1))] <- off_diag[mask]

    mask <- d_cmp > mp[off_diag]
    mp[c(rep(FALSE, diag - 1), mask)] <- d_cmp[mask]
    mpi[c(rep(FALSE, diag - 1), mask)] <- offset[mask]
  }
}

maximum_subsequence <- function(ts, threshold) {
  # """
  #   Finds the maximum subsequence length based on the threshold provided. Note
  #   that this threshold is domain specific requiring some knowledge about the
  #   underyling time series in question.
  #
  #   The subsequence length starts at 8 and iteratively doubles until the
  #   maximum correlation coefficient is no longer met.
  #
  #   Parameters
  #   ----------
  #   ts : array_like
  #       The time series to analyze.
  #   threshold : float
  #       The correlation coefficient used as the threshold. It should be between
  #       0 and 1.
  #
  #   Returns
  #   -------
  #   int :
  #       The maximum subsequence length based on the threshold provided.
  #   """
  #

  correlation_max <- Inf
  window_size <- 650
  max_window <- floor(length(ts) / 2)

  while (window_size <= max_window) {
    pmp <- mpx(ts, window_size)$mp
    mask <- !is.infinite(pmp)
    correlation_max <- max(pmp[mask])

    message("correlation_max: ", correlation_max, " window_size: ", window_size)

    if (correlation_max < threshold) {
      break
    }

    window_size <- window_size + 1
  }

  return(window_size)
}
