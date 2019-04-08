
`[.MultiMotif` <- function(x, ..., drop = FALSE) {
  x <- NextMethod(object = x)
  message("subMultiMotif")

  mp_size <- nrow(x$mp)
  offset <- attr(x, "offset")

  if (!is.null(offset) && offset > 0) {
    x$motif$motif_idx <- lapply(x$motif$motif_idx, function(y) {
      y - offset
    })

    x$motif$motif_neighbor <- lapply(x$motif$motif_neighbor, function(y) {
      y - offset
    })
  }

  x$motif$motif_idx <- lapply(x$motif$motif_idx, function(y) {
    y <- y[y <= mp_size & y > 0]

    if (length(y) < 2) {
      return(NULL)
    }
    y
  })

  # remove NULL's
  nulls <- sapply(x$motif$motif_idx, is.null)

  x$motif$motif_idx[nulls] <- NULL
  x$motif$motif_neighbor[nulls] <- NULL
  x$motif$motif_window[nulls] <- NULL
  x$motif$motif_neighbor <- lapply(x$motif$motif_neighbor, function(y) {
    y[y <= mp_size & y > 0]
  })

  x
}

`[.Motif` <- function(x, ..., drop = FALSE) {
  x <- NextMethod(object = x)
  message("subMotif")

  mp_size <- nrow(x$mp)
  offset <- attr(x, "offset")

  if (!is.null(offset) && offset > 0) {
    x$motif$motif_idx <- lapply(x$motif$motif_idx, function(y) {
      y - offset
    })

    x$motif$motif_neighbor <- lapply(x$motif$motif_neighbor, function(y) {
      y - offset
    })
  }

  x$motif$motif_idx <- lapply(x$motif$motif_idx, function(y) {
    y <- y[y <= mp_size & y > 0]

    if (length(y) < 2) {
      return(NULL)
    }
    y
  })

  # remove NULL's
  nulls <- sapply(x$motif$motif_idx, is.null)

  x$motif$motif_idx[nulls] <- NULL
  x$motif$motif_neighbor[nulls] <- NULL
  x$motif$motif_window[nulls] <- NULL
  x$motif$motif_neighbor <- lapply(x$motif$motif_neighbor, function(y) {
    y[y <= mp_size & y > 0]
  })

  x
}

`[.MatrixProfile` <- function(x, ..., drop = FALSE) {
  message("subMatrixProfile")
  # str(...)
  # y <- NextMethod("[")
  # y
  subset <- c(...)
  sub_size <- length(subset)

  if (is.numeric(subset)) {
    if (!all(diff(subset) == 1)) {
      stop("Indexes must be continuous and ascending.")
    }

    if (sub_size < (2 * max(x$w))) {
      stop("Subset must be larger than twice the window size: 2 * ", max(x$w))
    }

    if (!is.null(x$data)) {
      max_valid_idx <- nrow(x$data[[1]])
    } else {
      max_valid_idx <- nrow(x$mp) + min(x$w) - 1
    }

    if (max(subset) > max_valid_idx) {
      stop("Index is larger than data size.")
    } else if (max(subset) < max_valid_idx) {
      attr(x, "subset") <- TRUE
    }

    if (!is.null(x$data)) {
      # TODO: does query can be subset or just data?
      for (i in seq_along(x$data)) {
        x$data[[i]] <- x$data[[i]][subset, , drop = drop]
      }
    }

    mp_size <- sub_size - min(x$w) + 1
    mp_set <- subset[seq_len(mp_size)]
    offset <- subset[1] - 1

    x$mp <- x$mp[mp_set, , drop = drop]
    x$pi <- x$pi[mp_set, , drop = drop] - offset

    if (inherits(x, "Valmod")) {
      x$w <- x$w[mp_set, , drop = drop]
      x$mpnn <- x$mpnn[mp_set, , drop = drop]
      x$pinn <- x$pinn[mp_set, , drop = drop] - offset
      x$wnn <- x$wnn[mp_set, , drop = drop]
    }

    if (!is.null(attr(x, "join")) && !attr(x, "join") && !inherits(x, "Valmod")) {
      # rmp and lmp are not used in join-similarity
      x$rmp <- x$rmp[mp_set, , drop = drop]
      x$rpi <- x$rpi[mp_set, , drop = drop] - offset
      x$lmp <- x$lmp[mp_set, , drop = drop]
      x$lpi <- x$lpi[mp_set, , drop = drop] - offset
    }

    if (is.null(attr(x, "offset"))) {
      attr(x, "offset") <- offset
    } else {
      attr(x, "offset") <- attr(x, "offset") + offset
    }

    return(x)
  } else {
    getElement(x, args)
  }
}

`[.MultiMatrixProfile` <- function(x, ..., drop = FALSE) {
  message("subMultiMatrixProfile")

  subset <- c(...)
  sub_size <- length(subset)

  if (is.numeric(subset)) {
    if (!all(diff(subset) == 1)) {
      stop("Indexes must be continuous and ascending.")
    }

    if (sub_size < (2 * max(x$w))) {
      stop("Subset must be larger than twice the window size: 2 * ", max(x$w))
    }

    if (!is.null(x$data)) {
      max_valid_idx <- nrow(x$data[[1]])
    } else {
      max_valid_idx <- nrow(x$mp) + min(x$w) - 1
    }

    if (max(subset) > max_valid_idx) {
      stop("Index is larger than data size.")
    } else if (max(subset) < max_valid_idx) {
      attr(x, "subset") <- TRUE
    }

    if (!is.null(x$data)) {
      # TODO: does query can be subset or just data?
      for (i in seq_along(x$data)) {
        x$data[[i]] <- x$data[[i]][subset, , drop = drop]
      }
    }

    mp_size <- sub_size - min(x$w) + 1
    mp_set <- subset[seq_len(mp_size)]
    offset <- subset[1] - 1

    x$mp <- x$mp[mp_set, , drop = drop]
    x$pi <- x$pi[mp_set, , drop = drop] - offset

    if (!is.null(attr(x, "join")) && !attr(x, "join")) {
      # rmp and lmp are not used in join-similarity
      x$rmp <- x$rmp[mp_set, , drop = drop]
      x$rpi <- x$rpi[mp_set, , drop = drop] - offset
      x$lmp <- x$lmp[mp_set, , drop = drop]
      x$lpi <- x$lpi[mp_set, , drop = drop] - offset
    }

    if (is.null(attr(x, "offset"))) {
      attr(x, "offset") <- offset
    } else {
      attr(x, "offset") <- attr(x, "offset") + offset
    }

    return(x)
  } else {
    getElement(x, args)
  }
}

`[.SimpleMatrixProfile` <- function(x, ..., drop = FALSE) {
  message("subSimpleMatrixProfile")

  subset <- c(...)
  sub_size <- length(subset)

  if (is.numeric(subset)) {
    if (!all(diff(subset) == 1)) {
      stop("Indexes must be continuous and ascending.")
    }

    if (sub_size < (2 * max(x$w))) {
      stop("Subset must be larger than twice the window size: 2 * ", max(x$w))
    }

    if (!is.null(x$data)) {
      max_valid_idx <- nrow(x$data[[1]])
    } else {
      max_valid_idx <- nrow(x$mp) + min(x$w) - 1
    }

    if (max(subset) > max_valid_idx) {
      stop("Index is larger than data size.")
    } else if (max(subset) < max_valid_idx) {
      attr(x, "subset") <- TRUE
    }

    if (!is.null(x$data)) {
      # TODO: does query can be subset or just data?
      for (i in seq_along(x$data)) {
        x$data[[i]] <- x$data[[i]][subset, , drop = drop]
      }
    }

    mp_size <- sub_size - min(x$w) + 1
    mp_set <- subset[seq_len(mp_size)]
    offset <- subset[1] - 1

    x$mp <- x$mp[mp_set, , drop = drop]
    x$pi <- x$pi[mp_set, , drop = drop] - offset

    if (is.null(attr(x, "offset"))) {
      attr(x, "offset") <- offset
    } else {
      attr(x, "offset") <- attr(x, "offset") + offset
    }

    return(x)
  } else {
    getElement(x, args)
  }
}

tail.MatrixProfile <- function(x, n = 2 * max(x$w), ...) {
  data_size <- nrow(x$mp) + min(x$w) - 1

  if (n > 0) {
    st_idx <- data_size - n + 1
  } else {
    st_idx <- abs(n) + 1
  }

  message(paste0(st_idx, ":", data_size))
  return(x[st_idx:data_size])
}

tail.MultiMatrixProfile <- function(x, n = 2 * max(x$w), ...) {
  return(tail.MatrixProfile(x, n, ...))
}

tail.SimpleMatrixProfile <- function(x, n = 2 * max(x$w), ...) {
  return(tail.MatrixProfile(x, n, ...))
}

head.MatrixProfile <- function(x, n = 2 * max(x$w), ...) {
  data_size <- nrow(x$mp) + min(x$w) - 1

  if (n > 0) {
    ed_idx <- n
  } else {
    ed_idx <- data_size - abs(n)
  }

  message(paste0("1:", ed_idx))
  return(x[1:ed_idx])
}

head.MultiMatrixProfile <- function(x, n = 2 * max(x$w), ...) {
  return(head.MatrixProfile(x, n, ...))
}

head.SimpleMatrixProfile <- function(x, n = 2 * max(x$w), ...) {
  return(head.MatrixProfile(x, n, ...))
}
