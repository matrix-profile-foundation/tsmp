#' MDL Based MOTIF Discovery for Multidimensional Matrix Profile
#'
#' MDL Based MOTIF Discovery for Multidimensional Matrix Profile.
#'
#' Although this functions handles Multivariate Time Series, it can also be used to handle Univariate Time Series.
#'
#' @param data a `matrix` of `numeric`, where each colums is a time series. Accepts `vector` (see details), `list` and `data.frame` too.
#' @param window.size an `int` with the size of the sliding window.
#' @param matrix.profile multidimensional matrix profile (from [mstomp()] or [mstomp.par()]).
#' @param profile.index multidimensional profile index (from [mstomp()] or [mstomp.par()]).
#' @param n.bit an `int`. Number of bits for MDL discretization. (Default is `4`).
#' @param k an `int`. The number of MOTIFs to retrieve. `Inf` means all possible MOTIFs. (Default is `Inf`).
#'
#' @return Returns the `motif.idx` with the index of MOTIFs founded and `motif.dim`
#' with the spanned dimensions of respective MOTIF.
#' @export
#'
#' @family mstomp
#' @seealso [mstomp()], [mstomp.par()], [guide.search()]
#' @references 1. Yeh CM, Kavantzas N, Keogh E. Matrix Profile VI : Meaningful Multidimensional Motif Discovery.
#' @references 2. Zhu Y, Imamura M, Nikovski D, Keogh E. Matrix Profile VII: Time Series Chains: A New Primitive for Time Series Data Mining. Knowl Inf Syst. 2018 Jun 2;1–27.
#' @references Website: <https://sites.google.com/view/mstamp/>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#' \dontrun{
#' mp <- mstomp.par(toy_data$data, 30)
#' motifs <- unconstrain.search(toy_data$data, 30, mp$mp, mp$pi, 4, 2)
#' }
#'

unconstrain.search <- function(data, window.size, matrix.profile, profile.index, n.bit = 4, k = Inf) {

  ## transform data list into matrix
  if (is.list(data)) {
    data.len <- length(data[[1]])
    n.dim <- length(data)

    for (i in 1:n.dim) {
      len <- length(data[[i]])
      # Fix TS size with NaN
      if (len < data.len) {
        data[[i]] <- c(data[[i]], rep(NA, data.len - len))
      }
    }
    # transform data into matrix (each column is a TS)
    data <- sapply(data, cbind)
  } else if (is.matrix(data) || is.data.frame(data)) {
    if (is.data.frame(data)) {
      data <- as.matrix(data)
    } # just to be uniform
    if (ncol(data) > nrow(data))
      data <- t(data)
    data.len <- nrow(data)
    n.dim <- ncol(data)
  } else if (is.vector(data)) {
    data.len <- length(data)
    n.dim <- 1
    # transform data into 1-col matrix
    data <- as.matrix(data) # just to be uniform
  } else {
    stop("Unknown type of data. Must be: matrix, data.frame, vector or list")
  }

  exc.zone <- round(0.5 * window.size)
  tot.dim <- n.dim

  if (is.infinite(k)) {
    k <- dim(matrix.profile)[1]
  }

  motif.idx <- rep(0, k)
  motif.dim <- list()

  base.bit <- n.bit * tot.dim * window.size * 2
  found <- 0
  for (i in 1:k) {
    message(sprintf("Searching for motif (%d)", i))

    idx.1 <- apply(matrix.profile, 2, which.min) # sort by column
    val <- matrix.profile[cbind(idx.1, 1:ncol(matrix.profile))]

    if (any(is.infinite(val))) {
      motif.idx <- motif.idx[1:(k - 1)]
      motif.dim <- motif.dim[1:(k - 1)]
      break
    }

    bit.sz <- rep(0, tot.dim)
    idx.2 <- rep(0, tot.dim)

    dim <- list()

    for (j in 1:tot.dim) {
      idx.2[j] <- profile.index[idx.1[j], j]
      motif.1 <- data[idx.1[j]:(idx.1[j] + window.size - 1), ]
      motif.2 <- data[idx.2[j]:(idx.2[j] + window.size - 1), ]

      bits <- get.bit.save(motif.1, motif.2, j, n.bit)

      bit.sz[j] <- bits$bit.sz
      dim[[j]] <- bits$dim.id
    }

    min.idx <- which.min(bit.sz)
    best.bit <- bit.sz[min.idx]

    if (best.bit > (base.bit)) {
      if (i == 1)
        message("No motifs found")

      motif.idx <- motif.idx[1:(k - 1)]
      motif.dim <- motif.dim[1:(k - 1)]
      break
    } else
      found = found + 1

    motif.idx[i] <- idx.1[min.idx]
    motif.dim[[i]] <- dim[[min.idx]]

    st.idx <- max(1, motif.idx[i] - exc.zone)

    ed.idx <- min((dim(matrix.profile)[1]), motif.idx[i] + exc.zone)

    matrix.profile[st.idx:ed.idx, ] <- Inf
  }

  if (i != 1)
    message(sprintf("Found %d motifs", found))

  motif.dim <- motif.dim[motif.idx != 0]
  motif.idx <- motif.idx[motif.idx != 0]

  return(list(motif.idx = motif.idx, motif.dim = motif.dim))
}

get.bit.save <- function(motif.1, motif.2, n.dim, n.bit) {

  if (is.vector(motif.1))
    motif.1 <- as.matrix(motif.1)

  if (is.vector(motif.2))
    motif.2 <- as.matrix(motif.2)

  tot.dim <- dim(motif.1)[2]
  window.size <- dim(motif.1)[1]
  split.pt <- get.desc.split.pt(n.bit)
  disc.1 <- discretization(motif.1, split.pt)
  disc.2 <- discretization(motif.2, split.pt)

  dim.id <- sort(apply(abs(disc.1 - disc.2), 2, sum), index.return = TRUE)$ix
  dim.id <- dim.id[1:n.dim]
  motif.diff <- disc.1[, dim.id] - disc.2[, dim.id]
  n.val <- length(unique(as.vector(motif.diff)))

  bit.sz <- n.bit * (tot.dim * window.size * 2 - n.dim * window.size)
  bit.sz <- bit.sz + n.dim * window.size * log2(n.val) + n.val * n.bit

  return(list(bit.sz = bit.sz, dim.id = dim.id))
}

discretization <- function(motif, split.pt) {

  if (is.vector(motif))
    motif <- as.matrix(motif)

  dimmotif <- dim(motif)

  for (i in 1:dimmotif[2]) {
    motif[, i] <- (motif[, i] - mean(motif[, i])) / std(motif[, i])
  }

  disc <- matrix(0, dimmotif[1], dimmotif[2])

  for (i in 1:length(split.pt)) {
    disc[motif < split.pt[i] & disc == 0] <- i
  }

  disc[disc == 0] <- length(split.pt) + 1

  return(disc)
}


get.desc.split.pt <- function(n.bit) {
  split.pt <- stats::qnorm((1:((2^n.bit) - 1)) / (2^n.bit), 0, 1)
  return(split.pt)
}
