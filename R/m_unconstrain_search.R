#' Title
#'
#' @param data
#' @param sub_len
#' @param pro_mul
#' @param pro_idx
#' @param n_bit
#' @param k
#'
#' @return
#' @export
#'
#' @examples
unconstrain_search <- function(data, sub_len, pro_mul, pro_idx, n_bit, k) {

  ## transform data list into matrix
  if (is.list(data)) {
    data_len <- length(data[[1]])
    n_dim <- length(data)

    for (i in 1:n_dim) {
      len <- length(data[[i]])
      # Fix TS size with NaN
      if (len < data_len) {
        data[[i]] <- c(data[[i]], rep(NA, data_len - len))
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
    data_len <- nrow(data)
    n_dim <- ncol(data)
  } else if (is.vector(data)) {
    data_len <- length(data)
    n_dim <- 1
    # transform data into 1-col matrix
    data <- as.matrix(data) # just to be uniform
  } else {
    stop("Unknown type of data. Must be: matrix, data.frame, vector or list")
  }

  exc_zone <- round(0.5 * sub_len)
  tot_dim <- n_dim

  if (is.infinite(k)) {
    k <- dim(pro_mul)[1]
  }

  motif_idx <- rep(0, k)
  motif_dim <- list()

  base_bit <- n_bit * tot_dim * sub_len * 2
  found <- 0
  for (i in 1:k) {
    message(sprintf("Searching for motif (%d)", i))

    idx_1 <- apply(pro_mul, 2, which.min) # sort by column
    val <- pro_mul[cbind(idx_1, 1:ncol(pro_mul))]

    if (any(is.infinite(val))) {
      motif_idx <- motif_idx[1:(k - 1)]
      motif_dim <- motif_dim[1:(k - 1)]
      break
    }

    bit_sz <- rep(0, tot_dim)
    idx_2 <- rep(0, tot_dim)

    dim <- list()

    for (j in 1:tot_dim) {
      idx_2[j] <- pro_idx[idx_1[j], j]
      motif_1 <- data[idx_1[j]:(idx_1[j] + sub_len - 1), ]
      motif_2 <- data[idx_2[j]:(idx_2[j] + sub_len - 1), ]

      bits <- get_bit_save(motif_1, motif_2, j, n_bit)

      bit_sz[j] <- bits$bit_sz
      dim[[j]] <- bits$dim_id
    }

    min_idx <- which.min(bit_sz)
    best_bit <- bit_sz[min_idx]

    if (best_bit > (base_bit)) {
      if (i == 1)
        message("No motifs found")

      motif_idx <- motif_idx[1:(k - 1)]
      motif_dim <- motif_dim[1:(k - 1)]
      break
    } else
      found = found + 1

    motif_idx[i] <- idx_1[min_idx]
    motif_dim[[i]] <- dim[[min_idx]]

    st_idx <- max(1, motif_idx[i] - exc_zone)

    ed_idx <- min((dim(pro_mul)[1]), motif_idx[i] + exc_zone)

    pro_mul[st_idx:ed_idx, ] <- Inf
  }

  if (i != 1)
    message(sprintf("Found %d motifs", found))

  motif_dim <- motif_dim[motif_idx != 0]
  motif_idx <- motif_idx[motif_idx != 0]

  return(list(motif_idx = motif_idx, motif_dim = motif_dim))
}

get_bit_save <- function(motif_1, motif_2, n_dim, n_bit) {

  if (is.vector(motif_1))
    motif_1 <- as.matrix(motif_1)

  if (is.vector(motif_2))
    motif_2 <- as.matrix(motif_2)

  tot_dim <- dim(motif_1)[2]
  sub_len <- dim(motif_1)[1]
  split_pt <- get_desc_split_pt(n_bit)
  disc_1 <- discretization(motif_1, split_pt)
  disc_2 <- discretization(motif_2, split_pt)

  dim_id <- sort(apply(abs(disc_1 - disc_2), 2, sum), index.return = TRUE)$ix
  dim_id <- dim_id[1:n_dim]
  motif_diff <- disc_1[, dim_id] - disc_2[, dim_id]
  n_val <- length(unique(as.vector(motif_diff)))

  bit_sz <- n_bit * (tot_dim * sub_len * 2 - n_dim * sub_len)
  bit_sz <- bit_sz + n_dim * sub_len * log2(n_val) + n_val * n_bit

  return(list(bit_sz = bit_sz, dim_id = dim_id))
}

discretization <- function(motif, split_pt) {

  if (is.vector(motif))
    motif <- as.matrix(motif)

  dimmotif <- dim(motif)

  for (i in 1:dimmotif[2]) {
    motif[, i] <- (motif[, i] - mean(motif[, i])) / std(motif[, i])
  }

  disc <- matrix(0, dimmotif[1], dimmotif[2])

  for (i in 1:length(split_pt)) {
    disc[motif < split_pt[i] & disc == 0] <- i
  }

  disc[disc == 0] <- length(split_pt) + 1

  return(disc)
}


get_desc_split_pt <- function(n_bit) {
  split_pt <- qnorm((1:((2^n_bit) - 1)) / (2^n_bit), 0, 1)
  return(split_pt)
}
