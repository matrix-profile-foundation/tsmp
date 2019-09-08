#' Title
#'
#' @param Dataset
#' @param Range
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
#'   unlike PMP, which contains exact distances for all subsequences of all lengths.
#'
#'
#' @return
#' @export
#'
#' @examples
SKIMP <- function(Dataset, Range = seq.int(1, seq_end, jump)) {
  #   SKIMP generates the Pan Matrix Profile (PMP), a data structure which contains the Matrix
  #   Profile information for all subsequence lengths in a range. If not ran to completion, the
  #   instruction PMP(IDX) yields the approximation as a matrix
  #
  #   INPUTS Dataset - A nonempty, row or column vector, consisting of real values. Range   - A
  #   positive and increasing sequence of integers which enumerate the subsequence lengths used to
  #   create the Pan Matrix Profile
  #
  #   OUTPUTS PMP - Pan Matrix Profile approximation. Each row contains a single run of the Matrix
  #   Profile for some subsequence length. If ran to completion this variable stores the exact Pan
  #   Matrix Profile; otherwise, PMP(IDX) yields the approximation IDX - The permutation of the rows
  #   in the PMP which correspond to the PMP approximation
  #

  # assert(max(INPUTS[['Range']]) < length(Dataset) / 2, '[SKIMP] Error: The maximum subsequence
  # length explored must be less than half the size of the dataset.')
  #
  # INITIALIZE

  # Determine the order in which we will explore the subsequence lengths used to run the matrix profile
  SplitIDX <- BinarySplit(length(Range))

  # Stores the intermediate values of each Matrix Profile where PMP(i,:) = MatrixProfile with
  # subsequence length INPUTS[['Range']](SplitIDX(i))
  PMP <- matrix(Inf, nrow = length(Range), ncol = length(Dataset))

  # Stores row indices for image (Image is not a copy of the PMP but rather a permutation (with
  # repetition) of its rows
  IDX <- rep(-1, length(Range))

  ## BEGIN PROCESSING

  # Runs matrix profile on each subsequence length in INPUTS[['Range']](SplitIDX), stores the
  # results in PMP, and generates a frame for the animation
  for (i in 1:length(SplitIDX)) {

    # Get the current subsequence length to run the Matrix Profile on
    SubsequenceLength <- Range[SplitIDX[i]]

    # Run Matrix Profile
    # M <- mpx(Dataset, floor(SubsequenceLength / 4), SubsequenceLength)$mp
    M <- scrimp(Dataset, window_size = SubsequenceLength, pre_only = TRUE)$mp
    # to do ed
    PMP[SplitIDX[i], 1:length(M)] <- M

    # Propogate changes to IDX
    # TO-DO: Use a binary search to find last index, update all values in range. This MAY result in a time save
    j <- SplitIDX[i]
    while (j <= length(SplitIDX) && IDX[j] != j) {
      IDX[j] <- SplitIDX[i]
      j <- j + 1
    }

    test <- PMP[rev(IDX), ]
    depth <- 256
    test <- ceiling(test * depth) / depth
    test[test > 1] <- 1
    grid::grid.raster(test, height = 0.3, width = 1, interpolate = FALSE)

    message("i: ", i, "/", length(SplitIDX))
  }

  PMP <- PMP[rev(IDX), ]

  return(list(pmp = PMP, idx = IDX))
} # function SKIMP

#' Title
#'
#' @param n
#'
#' @return
#' @export
#'
#' @examples
BinarySplit <- function(n) {
  # BinarySplit(Length)
  # Preallocate memory to store the traversed indices and the remaining intervals
  #
  # Unit: milliseconds
  # expr      min       lq     mean   median       uq      max neval
  # BinarySplit(10000) 174.7164 180.7704 186.4403 182.1545 186.3419 234.9192   100
  #
  #
  IDX <- vector(mode = "numeric", length = n)
  Intervals <- list()

  IDX[1] <- 1 # We always begin by explore the first integer
  Intervals[[1]] <- c(2, n) # After exploring the first integer, we begin splitting the interval 2:n
  i <- 2

  while (length(Intervals) > 0) {
    lb <- Intervals[[1]][1]
    ub <- Intervals[[1]][2]
    mid <- floor((lb + ub) / 2)
    Intervals[[1]] <- NULL

    IDX[i] <- mid
    i <- i + 1

    if (lb == ub) {
      next
    } else {
      lr <- split(lb, ub, mid)
      if (!is.null(lr$L)) {
        Intervals[[length(Intervals) + 1]] <- lr$L
      }
      if (!is.null(lr$R)) {
        Intervals[[length(Intervals) + 1]] <- lr$R
      }
    }
  }
  return(IDX)
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
    L <- NULL
    R <- c(m + 1, ub)
  } else if (ub == m) {
    L <- c(lb, m - 1)
    R <- NULL
  } else {
    L <- c(lb, m - 1)
    R <- c(m + 1, ub)
  }

  return(list(L = L, R = R))
}



## Functions from various sources brought together into a single document with various authors

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
mpx <- function(a, minlag, w) {

  # matrix profile using cross correlation,

  # depends on files sum2s, musigtest, dot2s
  tictac <- Sys.time()
  n <- length(a)

  msd <- fast_muinvn(a, w)
  mu <- msd$mu
  sig <- msd$sig

  # differentials have 0 as their first entry. This simplifies index
  # calculations slightly and allows us to avoid special "first line"
  # handling.
  df <- c(0, (1 / 2) * (a[(1 + w):n] - a[1:(n - w)]))
  dg <- c(0, (a[(1 + w):n] - mu[2:(n - w + 1)]) + (a[1:(n - w)] - mu[1:(n - w)]))
  diagmax <- length(a) - w + 1
  mp <- kronecker(matrix(1, n - w + 1, 1), -1)
  mpi <- rep(NA, n - w + 1)

  for (diag in (minlag + 1):diagmax) {
    c <- sum((a[diag:(diag + w - 1)] - mu[diag]) * (a[1:w] - mu[1]))
    for (offset in 1:(n - w - diag + 2)) {
      c <- c + df[offset] * dg[offset + diag - 1] + df[offset + diag - 1] * dg[offset]
      c_cmp <- c * sig[offset] * sig[offset + diag - 1]
      if (c_cmp > mp[offset]) {
        mp[offset] <- c_cmp
        mpi[offset] <- offset + diag - 1
      }
      if (c_cmp > mp[offset + diag - 1]) {
        mp[offset + diag - 1] <- c_cmp
        mpi[offset + diag - 1] <- offset
      }
    }
  }
  # to do ed
  mp <- sqrt(2 * w * (1 - mp))
  tictac <- Sys.time() - tictac

  message(sprintf("mpx finished in %.2f %s", tictac, units(tictac)))

  return(list(mp = mp, mpi = mpi))
}
