#---- Subset Chains ----

#' @export
#' @keywords internal
#' @noRd
#'

`[.Chain` <- function(x, ..., drop = FALSE) {
  x <- NextMethod(object = x)
  attr(x, "subsetting") <- "Chains"

  mp_size <- nrow(x$mp)
  offset <- attr(x, "offset")

  if (is.null(offset)) {
    offset <- 0
  }

  if (offset > 0) {
    x$chain$chains <- lapply(x$chain$chains, function(y) {
      y - offset
    })

    x$chain$best <- x$chain$best - offset
  }

  x$chain$chains <- lapply(x$chain$chains, function(y) {
    y <- y[y <= mp_size & y > 0]

    if (length(y) < 3) {
      return(NULL)
    }
    y
  })

  x$chain$best <- x$chain$best[x$chain$best <= mp_size & x$chain$best > 0]

  # remove NULL's
  nulls <- sapply(x$chain$chains, is.null)
  x$chain$chains[nulls] <- NULL

  attr(x, "subsetting") <- NULL

  x
}

#---- Subset Salient ----

#' @export
#' @keywords internal
#' @noRd
#'

`[.Salient` <- function(x, ..., drop = FALSE) {
  x <- NextMethod(object = x)
  attr(x, "subsetting") <- "Salient"

  offset <- attr(x, "offset")

  if (is.null(offset)) {
    offset <- 0
  }

  x$salient$indexes <- x$salient$indexes - offset
  idxs <- (x$salient$indexes > 0 & x$salient$indexes <= nrow(x$mp))
  x$salient$indexes <- as.matrix(x$salient$indexes[idxs])
  x$salient$idx_bit_size <- as.matrix(x$salient$idx_bit_size[idxs])

  attr(x, "subsetting") <- NULL

  x
}

#---- Subset Annotations ----

#' @export
#' @keywords internal
#' @noRd
#'

`[.AnnotationVector` <- function(x, ..., drop = FALSE) {
  x <- NextMethod(object = x)
  attr(x, "subsetting") <- "AnnotationVector"

  offset <- attr(x, "offset")

  if (is.null(offset)) {
    offset <- 0
  }

  st_idx <- 1 + offset
  ed_idx <- st_idx + nrow(x$mp) - 1
  x$av <- x$av[st_idx:ed_idx]

  attr(x, "subsetting") <- NULL

  x
}

#---- Subset Fluss ----

#' @export
#' @keywords internal
#' @noRd
#'

`[.Fluss` <- function(x, ..., drop = FALSE) {
  x <- NextMethod(object = x)
  attr(x, "subsetting") <- "Fluss"

  x <- fluss_extract(x, length(x$fluss))

  attr(x, "subsetting") <- NULL

  x
}

#' @export
#' @keywords internal
#' @noRd
#'

`[.Floss` <- function(x, ..., drop = FALSE) {
  x <- NextMethod(object = x)
  attr(x, "subsetting") <- "Floss"

  x <- floss_extract(x, length(x$floss))

  attr(x, "subsetting") <- NULL

  x
}

#' @export
#' @keywords internal
#' @noRd
#'

`[.ArcCount` <- function(x, ..., drop = FALSE) {
  x <- NextMethod(object = x)
  attr(x, "subsetting") <- "ArcCount"

  x <- fluss_cac(x)

  attr(x, "subsetting") <- NULL

  x
}

#---- Subset Discord ----

#' @export
#' @keywords internal
#' @noRd
#'

`[.Discord` <- function(x, ..., drop = FALSE) {
  x <- NextMethod(object = x)
  attr(x, "subsetting") <- "Discord"

  mp_size <- nrow(x$mp)
  offset <- attr(x, "offset")

  if (!is.null(offset) && offset > 0) {
    x$discord$discord_idx <- lapply(x$discord$discord_idx, function(y) {
      y - offset
    })

    x$discord$discord_neighbor <- lapply(x$discord$discord_neighbor, function(y) {
      y - offset
    })
  }

  x$discord$discord_idx <- lapply(x$discord$discord_idx, function(y) {
    y <- y[y <= mp_size & y > 0]

    if (length(y) < 1) {
      return(NULL)
    }
    y
  })

  # remove NULL's
  nulls <- sapply(x$discord$discord_idx, is.null)

  x$discord$discord_idx[nulls] <- NULL
  x$discord$discord_neighbor[nulls] <- NULL
  x$discord$discord_window[nulls] <- NULL
  x$discord$discord_neighbor <- lapply(x$discord$discord_neighbor, function(y) {
    y[y <= mp_size & y > 0]
  })

  attr(x, "subsetting") <- NULL

  x
}

#---- Subset Motifs ----

#' @export
#' @keywords internal
#' @noRd
#'

`[.MultiMotif` <- function(x, ..., drop = FALSE) {
  x <- NextMethod(object = x)
  attr(x, "subsetting") <- "MultiMotif"

  mp_size <- nrow(x$mp)
  offset <- attr(x, "offset")

  # adjust offset
  if (!is.null(offset) && offset > 0) {
    x$motif$motif_idx <- lapply(x$motif$motif_idx, function(y) {
      y - offset
    })

    x$motif$motif_neighbor <- lapply(x$motif$motif_neighbor, function(y) {
      y - offset
    })
  }

  # remove from list
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
  x$motif$motif_dim[nulls] <- NULL
  x$motif$motif_neighbor <- lapply(x$motif$motif_neighbor, function(y) {
    y[y <= mp_size & y > 0]
  })

  attr(x, "subsetting") <- NULL
  x
}

#' @export
#' @keywords internal
#' @noRd
#'

`[.Motif` <- function(x, ..., drop = FALSE) {
  x <- NextMethod(object = x)
  attr(x, "subsetting") <- "Motif"

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

  attr(x, "subsetting") <- NULL

  x
}

#---- Subset Matrices ----

#' @export
#' @keywords internal
#' @noRd
#'

`[.PMP` <- function(x, ..., drop = FALSE) {
  stop("Subsetting PMP is not implemented yet")
}

#' @export
#' @keywords internal
#' @noRd
#'

`[.MatrixProfile` <- function(x, ..., drop = FALSE) {
  # str(...)
  # y <- NextMethod("[")
  # y
  subset <- c(...)
  sub_size <- length(subset)

  if (is.numeric(subset)) {
    attr(x, "subsetting") <- "MatrixProfile"

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
      new_data <- attr(x, "new_data")

      if (!is.null(new_data)) {
        removed <- nrow(x$mp) - max(subset)

        if (new_data <= removed) {
          attr(x, "new_data") <- 0
        } else {
          attr(x, "new_data") <- new_data - removed
        }
      }
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

    attr(x, "subsetting") <- NULL
    return(x)
  } else {
    getElement(x, args)
  }
}

#' @export
#' @keywords internal
#' @noRd
#'

`[.MultiMatrixProfile` <- function(x, ..., drop = FALSE) {
  subset <- c(...)
  sub_size <- length(subset)

  if (is.numeric(subset)) {
    attr(x, "subsetting") <- "MultiMatrixProfile"

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

    attr(x, "subsetting") <- NULL

    return(x)
  } else {
    getElement(x, args)
  }
}

#' @export
#' @keywords internal
#' @noRd
#'

`[.SimpleMatrixProfile` <- function(x, ..., drop = FALSE) {
  subset <- c(...)
  sub_size <- length(subset)

  if (is.numeric(subset)) {
    attr(x, "subsetting") <- "SimpleMatrixProfile"

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

    attr(x, "subsetting") <- NULL

    return(x)
  } else {
    getElement(x, args)
  }
}

#---- Subset MPdist ----

`[.MPdistProfile` <- function(x, ..., drop = FALSE) {
  subset <- c(...)
  sub_size <- length(subset)

  if (is.numeric(subset)) {
    attr(x, "subsetting") <- "MPdistProfile"

    if (!all(diff(subset) == 1)) {
      stop("Indexes must be continuous and ascending.")
    }

    attr <- attr(x, "origin")

    if (sub_size < attr$query_size) {
      stop("Subset must be at least the size of query data", attr$query_size)
    }

    if (!is.null(x$data)) {
      max_valid_idx <- nrow(x$data[[1]])
    } else {
      max_valid_idx <- nrow(x$mpdist) - attr$query_size + 1
    }

    if (max(subset) > max_valid_idx) {
      stop("Index is larger than data size.")
    } else if (max(subset) < max_valid_idx) {
      attr(x, "subset") <- TRUE

      # TODO: Handle new data
      # new_data <- attr(x, "new_data")
      #
      # if (!is.null(new_data)) {
      #   removed <- nrow(x$mp) - max(subset)
      #
      #   if (new_data <= removed) {
      #     attr(x, "new_data") <- 0
      #   } else {
      #     attr(x, "new_data") <- new_data - removed
      #   }
      # }
    }

    if (!is.null(x$data)) {
      x$data[[1]] <- x$data[[1]][subset, , drop = drop]
    }

    mp_size <- sub_size - attr$query_size + 1
    mp_set <- subset[seq_len(mp_size)]
    offset <- subset[1] - 1

    x$mpdist <- x$mpdist[mp_set, , drop = drop]

    if (is.null(attr(x, "offset"))) {
      attr(x, "offset") <- offset
    } else {
      attr(x, "offset") <- attr(x, "offset") + offset
    }

    attr(x, "subsetting") <- NULL
    return(x)
  } else {
    getElement(x, args)
  }
}

`[.Snippet` <- function(x, ..., drop = FALSE) {
  attr(x, "subsetting") <- "Snippet"

  message("Not implemented yet.")

  attr(x, "subsetting") <- NULL

  x
}

#---- Tails ----

tail.MatrixProfile <- function(x, n = 2 * max(x$w), ...) {
  data_size <- nrow(x$mp) + min(x$w) - 1

  if (n > 0) {
    st_idx <- data_size - n + 1
  } else {
    st_idx <- abs(n) + 1
  }
  return(x[st_idx:data_size])
}

tail.MultiMatrixProfile <- function(x, n = 2 * max(x$w), ...) {
  return(tail.MatrixProfile(x, n, ...))
}

tail.SimpleMatrixProfile <- function(x, n = 2 * max(x$w), ...) {
  return(tail.MatrixProfile(x, n, ...))
}

#---- Heads ----

head.MatrixProfile <- function(x, n = 2 * max(x$w), ...) {
  data_size <- nrow(x$mp) + min(x$w) - 1

  if (n > 0) {
    ed_idx <- n
  } else {
    ed_idx <- data_size - abs(n)
  }
  return(x[1:ed_idx])
}


head.MultiMatrixProfile <- function(x, n = 2 * max(x$w), ...) {
  return(head.MatrixProfile(x, n, ...))
}


head.SimpleMatrixProfile <- function(x, n = 2 * max(x$w), ...) {
  return(head.MatrixProfile(x, n, ...))
}
