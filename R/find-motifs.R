#' Search for Motifs
#'
#' @param .mp a `MatrixProfile` or `MultiMatrixProfile` object.
#' @param \dots further arguments to be passed to class specific function.
#'
#' @name find_motif
#' @export

find_motif <- function(.mp, ...) {
  UseMethod("find_motif", .mp)
}

#' @param data the data used to build the Matrix Profile, if not embedded.
#' @param n_motifs an `int`. Number of motifs to find. (Default is `3`).
#' @param n_neighbors an `int`. Number of neighbors to find. (Default is `10`).
#' @param radius an `int`. Set a threshold to exclude matching neighbors with distance > current
#' motif distance * `radius`. (Default is `3`).
#' @param exclusion_zone if a `number` will be used instead of embedded value. (Default is `NULL`).
#'
#' @name find_motif
#' @export
#' @return For class `MatrixProfile`, returns the input `.mp` object with a new name `motif`. It contains: `motif_idx`, a `list`
#' of motif pairs found and `motif_neighbor` a `list` with respective motif's neighbors.
#' @examples
#' # Single dimension data
#' w <- 50
#' data <- mp_gait_data
#' mp <- tsmp(data, window_size = w, exclusion_zone = 1 / 4, verbose = 0)
#' mp <- find_motif(mp)
find_motif.MatrixProfile <- function(.mp, data, n_motifs = 3, n_neighbors = 10, radius = 3, exclusion_zone = NULL, ...) {
  if (!("MatrixProfile" %in% class(.mp))) {
    stop("First argument must be an object of class `MatrixProfile`.")
  }

  if (inherits(.mp, "Valmod")) {
    valmod <- TRUE
  } else {
    valmod <- FALSE
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
    data_len <- nrow(data)
    data_dim <- ncol(data)
  } else if (is.list(data)) {
    data_len <- length(data[[1]])
    data_dim <- length(data)

    for (i in 1:data_dim) {
      len <- length(data[[i]])
      # Fix TS size with NaN
      if (len < data_len) {
        data[[i]] <- c(data[[i]], rep(NA, data_len - len))
      }
    }
    # transform data into matrix (each column is a TS)
    data <- sapply(data, cbind)
  } else if (is.vector(data)) {
    data_len <- length(data)
    data_dim <- 1
    # transform data into 1-col matrix
    data <- as.matrix(data) # just to be uniform
  } else {
    stop("`data` must be `matrix`, `data.frame`, `vector` or `list`.")
  }


  matrix_profile <- .mp # keep mp intact
  matrix_profile_size <- length(matrix_profile$mp)
  motif_idxs <- list(motifs = list(NULL), neighbors = list(NULL), windows = list(NULL))

  if (is.null(exclusion_zone)) {
    exclusion_zone <- .mp$ez
  }

  exclusion_zone <- round(.mp$w * exclusion_zone + vars()$eps)

  nn <- NULL

  for (i in seq_len(n_motifs)) {
    idxs <- min_mp_idx(matrix_profile)

    if (is.na(idxs[1])) {
      break
    }

    min_idx <- idxs[1]
    motif_distance <- matrix_profile$mp[min_idx]
    motif_idxs[[1]][[i]] <- sort(idxs)
    motif_idx <- motif_idxs[[1L]][[i]][1]

    if (valmod) {
      # precompute for each window size in valmod
      nn <- NULL
      window <- .mp$w[min_idx]
      e_zone <- exclusion_zone[min_idx]
    } else {
      window <- .mp$w
      e_zone <- exclusion_zone
    }

    # query using the motif to find its neighbors
    nn <- dist_profile(data, data, nn, window_size = window, index = min_idx)

    distance_profile <- nn$distance_profile

    if (valmod) {
      distance_profile <- distance_profile * sqrt(1.0 / window)
    }

    distance_profile[distance_profile > (motif_distance * radius)^2] <- Inf
    motif_zone_start <- pmax(1, motif_idx - e_zone)
    motif_zone_end <- pmin(matrix_profile_size, motif_idx + e_zone)
    distance_profile[motif_zone_start:motif_zone_end] <- Inf
    motif_idx <- motif_idxs[[1]][[i]][2]
    motif_zone_start <- pmax(1, motif_idx - e_zone)
    motif_zone_end <- pmin(matrix_profile_size, motif_idx + e_zone)
    distance_profile[motif_zone_start:motif_zone_end] <- Inf
    st <- sort(distance_profile, index.return = TRUE)
    distance_order <- st$x
    distance_idx_order <- st$ix

    motif_neighbor <- vector(mode = "numeric")

    for (j in seq_len(n_neighbors)) {
      if (is.infinite(distance_order[1]) || length(distance_order) < j) {
        break
      }
      motif_neighbor[j] <- distance_idx_order[1]
      distance_order <- distance_order[2:length(distance_order)]
      distance_idx_order <- distance_idx_order[2:length(distance_idx_order)]
      distance_order <- distance_order[!(abs(distance_idx_order - motif_neighbor[j]) < e_zone)]
      distance_idx_order <- distance_idx_order[!(abs(distance_idx_order - motif_neighbor[j]) < e_zone)]
    }

    motif_neighbor <- motif_neighbor[motif_neighbor != 0]
    motif_idxs[[2]][[i]] <- motif_neighbor
    motif_idxs[[3]][[i]] <- window

    remove_idx <- c(motif_idxs[[1]][[i]], motif_idxs[[2]][[i]])

    for (j in seq_len(length(remove_idx))) {
      remove_zone_start <- max(1, remove_idx[j] - e_zone)
      remove_zone_end <- min(matrix_profile_size, remove_idx[j] + e_zone)
      matrix_profile$mp[remove_zone_start:remove_zone_end] <- Inf
    }
  }

  if (is.null(motif_idxs[[1]][[1]])) {
    message("No valid motif found.")
    .mp <- remove_class(.mp, "Motif")
    return(.mp)
  }

  .mp$motif <- list(motif_idx = motif_idxs[[1]], motif_neighbor = motif_idxs[[2]], motif_window = motif_idxs[[3]])
  class(.mp) <- update_class(class(.mp), "Motif")
  return(.mp)
}

#' @param mode a `string`. Guided or Unconstrained search. Allow partial match. (Default is `guided`).
#' @param n_bit an `int`. Bit size for discretization. Ignored on Guided search. (Default is `4`).
#' @param n_dim an `int`. Number of dimensions to use on Guided search instead of embedded value. (Default is `NULL`).
#'
#' @return For class `MultiMatrixProfile`, returns the input `.mp` object with a new name `motif`. It contains: `motif_idx`, a `vector`
#' of motifs found and `motif_dim` a `list` the dimensions where the motifs were found
#'
#' @name find_motif
#' @export
#' @examples
#'
#' # Multidimension data
#' w <- mp_toy_data$sub_len
#' data <- mp_toy_data$data[1:200, ]
#' mp <- tsmp(data, window_size = w, mode = "mstomp", verbose = 0)
#' mp <- find_motif(mp)
find_motif.MultiMatrixProfile <- function(.mp, data, n_motifs = 3, mode = c("guided", "unconstrained"),
                                          n_bit = 4, exclusion_zone = NULL, n_dim = NULL, ...) {
  if (!("MultiMatrixProfile" %in% class(.mp))) {
    stop("First argument must be an object of class `MultiMatrixProfile`.")
  }

  if (missing(data) && !is.null(.mp$data)) {
    data <- .mp$data[[1]]
  }

  algo <- match.arg(mode)

  # transform data list into matrix
  if (is.matrix(data) || is.data.frame(data)) {
    if (is.data.frame(data)) {
      data <- as.matrix(data)
    } # just to be uniform
    if (ncol(data) > nrow(data)) {
      data <- t(data)
    }
    data_len <- nrow(data)
    data_dim <- ncol(data)
  } else if (is.list(data)) {
    data_len <- length(data[[1]])
    data_dim <- length(data)

    for (i in 1:data_dim) {
      len <- length(data[[i]])
      # Fix TS size with NaN
      if (len < data_len) {
        data[[i]] <- c(data[[i]], rep(NA, data_len - len))
      }
    }
    # transform data into matrix (each column is a TS)
    data <- sapply(data, cbind)
  } else if (is.vector(data)) {
    data_len <- length(data)
    data_dim <- 1
    # transform data into 1-col matrix
    data <- as.matrix(data) # just to be uniform
  } else {
    stop("`data` must be `matrix`, `data.frame`, `vector` or `list`.")
  }

  # Guided Search ------------------------------------------------------------------------
  if (algo == "guided") {
    if (is.null(n_dim)) {
      if (.mp$n_dim != data_dim) {
        warning("Warning: `data` dimensions are different from matrix profile.")
      }
      n_dim <- .mp$n_dim
    }

    motif_idx <- sort(min_mp_idx(.mp, n_dim))

    if (!is.na(motif_idx[1])) {
      motif_1 <- as.matrix(data[motif_idx[1]:(motif_idx[1] + .mp$w - 1), ]) # as.matrix(): hack for vectors
      motif_2 <- as.matrix(data[motif_idx[2]:(motif_idx[2] + .mp$w - 1), ]) # as.matrix(): hack for vectors

      motif_dim <- sort(apply(abs(motif_1 - motif_2), 2, sum), index.return = TRUE)$ix
      motif_dim <- sort(motif_dim[1:n_dim])
      motif_dim <- list(motif_dim)
      motif_idx <- list(motif_idx)
    } else {
      motif_idx <- NULL
      motif_dim <- NULL
    }

    .mp$motif <- list(motif_idx = motif_idx, motif_neighbor = NULL, motif_dim = motif_dim)
    class(.mp) <- update_class(class(.mp), "MultiMotif")
    return(.mp)
  } else {
    # Unguided Search -------------------------------------------------------------------
    if (n_bit < 2) {
      stop("`nbit` must be at least `2`.")
    }

    if (is.null(exclusion_zone)) {
      exclusion_zone <- .mp$ez
    }
    exclusion_zone <- round(exclusion_zone * .mp$w + vars()$eps)
    matrix_profile <- .mp # keep mp intact

    if (.mp$n_dim != data_dim) {
      warning("Warning: `data` dimensions are different from matrix profile.")
    }

    tot_dim <- .mp$n_dim

    if (is.infinite(n_motifs)) {
      n_motifs <- dim(matrix_profile$mp)[1]
    }

    motif_idx <- list()
    motif_dim <- list()

    base_bit <- n_bit * tot_dim * .mp$w * 2
    found <- 0
    for (i in seq_len(n_motifs)) {
      message(sprintf("Searching for motif (%d).", i))

      idxs <- min_mp_idx(matrix_profile)
      if (is.na(idxs[1])) {
        motif_dim <- motif_dim[seq_len(i - 1)]
        break
      }

      val <- matrix_profile$mp[cbind(idxs[, 1], seq_len(ncol(matrix_profile$mp)))]

      if (any(is.infinite(val))) {
        motif_dim <- motif_dim[seq_len(i - 1)]
        break
      }

      bit_sz <- rep(Inf, tot_dim)
      dim <- list()

      for (j in seq_len(tot_dim)) {
        motif_1 <- data[idxs[j, 1]:(idxs[j, 1] + .mp$w - 1), ]
        motif_2 <- data[idxs[j, 2]:(idxs[j, 2] + .mp$w - 1), ]

        bits <- get_bit_save(motif_1, motif_2, j, n_bit)

        bit_sz[j] <- bits$bit_sz
        dim[[j]] <- bits$dim_id
      }

      min_idx <- which.min(bit_sz)
      best_bit <- bit_sz[min_idx]

      if (best_bit > (base_bit)) {
        if (i == 1) {
          message("No valid motif found.")
          .mp <- remove_class(.mp, "MultiMotif")
          return(.mp)
        }
        motif_idx <- motif_idx[seq_len(i - 1)]
        motif_dim <- motif_dim[seq_len(i - 1)]
        break
      } else {
        found <- found + 1
      }

      motif_idx[[found]] <- sort(c(idxs[min_idx, 1], idxs[min_idx, 2]))
      motif_dim[[i]] <- sort(dim[[min_idx]])

      st_idx <- max(1, idxs[min_idx, 1] - exclusion_zone)
      ed_idx <- min((dim(matrix_profile$mp)[1]), idxs[min_idx, 1] + exclusion_zone)

      matrix_profile$mp[st_idx:ed_idx, ] <- Inf

      st_idx <- max(1, idxs[min_idx, 2] - exclusion_zone)
      ed_idx <- min((dim(matrix_profile$mp)[1]), idxs[min_idx, 2] + exclusion_zone)

      matrix_profile$mp[st_idx:ed_idx, ] <- Inf
    }

    motif_dim # <- motif_dim[motif_idx != 0]

    if (length(motif_idx) > 0) {
      message(sprintf("Found %d motifs.", found))
    }

    .mp$motif <- list(motif_idx = motif_idx, motif_neighbor = NULL, motif_dim = motif_dim)
    class(.mp) <- update_class(class(.mp), "MultiMotif")
    return(.mp)
  }
}

#' @param data the data used to build the Matrix Profile, if not embedded.
#' @param n_motifs an `int`. Number of motifs to find. (Default is `3`).
#' @param n_neighbors an `int`. Number of neighbors to find. (Default is `10`).
#' @param radius an `int`. Set a threshold to exclude matching neighbors with distance > current
#' motif distance * `radius`. (Default is `3`).
#' @param exclusion_zone if a `number` will be used instead of embedded value. (Default is `NULL`).
#'
#' @name find_motif
#' @export
#' @return For class `PMP`, returns the input `.mp` object with a new name `motif`. It contains: `motif_idx`, a `list`
#' of motif pairs found and `motif_neighbor` a `list` with respective motif's neighbors.
#' @examples
#' pan <- tsmp(mp_gait_data, window_size = 20:30, mode = "pmp")
#' mp <- find_motif(pan)
find_motif.PMP <- function(.mp, data, n_motifs = 3, n_neighbors = 10, radius = 3, exclusion_zone = NULL, ...) {
  if (!("PMP" %in% class(.mp))) {
    stop("First argument must be an object of class `PMP`.")
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
    data_len <- nrow(data)
    data_dim <- ncol(data)
  } else if (is.list(data)) {
    data_len <- length(data[[1]])
    data_dim <- length(data)

    for (i in 1:data_dim) {
      len <- length(data[[i]])
      # Fix TS size with NaN
      if (len < data_len) {
        data[[i]] <- c(data[[i]], rep(NA, data_len - len))
      }
    }
    # transform data into matrix (each column is a TS)
    data <- sapply(data, cbind)
  } else if (is.vector(data)) {
    data_len <- length(data)
    data_dim <- 1
    # transform data into 1-col matrix
    data <- as.matrix(data) # just to be uniform
  } else {
    stop("`data` must be `matrix`, `data.frame`, `vector` or `list`.")
  }

  # TODO: for now, computes only for the first window:
  row <- 1

  matrix_profile <- .mp # keep mp intact
  matrix_profile_size <- length(matrix_profile$pmp[[row]])
  motif_idxs <- list(motifs = list(NULL), neighbors = list(NULL), windows = list(NULL))

  if (is.null(exclusion_zone)) {
    exclusion_zone <- .mp$ez
  }

  exclusion_zone <- round(.mp$w[row] * exclusion_zone + vars()$eps)

  nn <- NULL

  for (i in seq_len(n_motifs)) {
    # idxs <- min_mp_idx(matrix_profile) # todo for PMP
    idxs <- which.min(matrix_profile$pmp[[row]])
    idxs <- c(idxs, matrix_profile$pmpi[[row]][idxs])


    if (is.na(idxs[1])) {
      break
    }

    min_idx <- idxs[1]
    motif_distance <- matrix_profile$pmp[[row]][min_idx]
    motif_idxs[[1]][[i]] <- sort(idxs)
    motif_idx <- motif_idxs[[1L]][[i]][1]

    window <- .mp$w[1]
    e_zone <- exclusion_zone


    # query using the motif to find its neighbors
    nn <- dist_profile(data, data, nn, window_size = window, index = min_idx)

    distance_profile <- nn$distance_profile

    distance_profile[distance_profile > (motif_distance * radius)^2] <- Inf
    motif_zone_start <- pmax(1, motif_idx - e_zone)
    motif_zone_end <- pmin(matrix_profile_size, motif_idx + e_zone)
    distance_profile[motif_zone_start:motif_zone_end] <- Inf
    motif_idx <- motif_idxs[[1]][[i]][2]
    motif_zone_start <- pmax(1, motif_idx - e_zone)
    motif_zone_end <- pmin(matrix_profile_size, motif_idx + e_zone)
    distance_profile[motif_zone_start:motif_zone_end] <- Inf
    st <- sort(distance_profile, index.return = TRUE)
    distance_order <- st$x
    distance_idx_order <- st$ix

    motif_neighbor <- vector(mode = "numeric")

    for (j in seq_len(n_neighbors)) {
      if (is.infinite(distance_order[1]) || length(distance_order) < j) {
        break
      }
      motif_neighbor[j] <- distance_idx_order[1]
      distance_order <- distance_order[2:length(distance_order)]
      distance_idx_order <- distance_idx_order[2:length(distance_idx_order)]
      distance_order <- distance_order[!(abs(distance_idx_order - motif_neighbor[j]) < e_zone)]
      distance_idx_order <- distance_idx_order[!(abs(distance_idx_order - motif_neighbor[j]) < e_zone)]
    }

    motif_neighbor <- motif_neighbor[motif_neighbor != 0]
    motif_idxs[[2]][[i]] <- motif_neighbor
    motif_idxs[[3]][[i]] <- window

    remove_idx <- c(motif_idxs[[1]][[i]], motif_idxs[[2]][[i]])

    for (j in seq_len(length(remove_idx))) {
      remove_zone_start <- max(1, remove_idx[j] - e_zone)
      remove_zone_end <- min(matrix_profile_size, remove_idx[j] + e_zone)
      matrix_profile$pmp[[row]][remove_zone_start:remove_zone_end] <- Inf
    }
  }

  if (is.null(motif_idxs[[1]][[1]])) {
    message("No valid motif found.")
    .mp <- remove_class(.mp, "Motif")
    return(.mp)
  }

  .mp$motif <- list(motif_idx = motif_idxs[[1]], motif_neighbor = motif_idxs[[2]], motif_window = motif_idxs[[3]])
  class(.mp) <- update_class(class(.mp), "Motif")
  return(.mp)
}
