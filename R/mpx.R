
#' Fast implementation of MP and MPI for internal purposes, without FFT
#'
#' @param a data
#' @param w window size
#' @param idx calc and return indexes?
#' @param rcpp uses rcpp?
#'
#' @return Returns MP and MPI
#' @keywords internal
#' @noRd

mpx <- function(a, w, idx = FALSE, dist = c("euclidean", "pearson"), rcpp = TRUE) {
  minlag <- floor(w / 4)
  dist <- match.arg(dist)

  if (dist == "euclidean") {
    dist <- TRUE
  } else {
    dist <- FALSE
  }

  if (rcpp) {
    return(mpx_rcpp(
      as.numeric(a),
      as.integer(w),
      as.integer(minlag),
      as.logical(idx),
      as.logical(dist)
    ))
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

  # anytime must return the result always
  on.exit(return({
    # to do ed
    mp[mp > 1] <- 1.0

    if (dist) {
      mp <- sqrt(2 * w * (1 - mp))
    }
    list(mp = mp, pi = mpi)
  }), TRUE)

  bb <- (a[1:w] - mu[1])

  for (diag in seq_diag) {
    c <- (a[diag:(diag + w - 1)] - mu[diag]) %*% bb

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
