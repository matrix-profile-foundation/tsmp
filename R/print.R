#' Prints a Valmod Matrix Profile
#'
#' @param x a TSMP object of class `Valmod`.
#' @param \dots additional arguments ignored.
#' @export
#' @keywords internal
#' @noRd

print.Valmod <- function(x, ...) {
  cat("Valmod Matrix Profile\n")
  cat("---------------------\n")

  cat("Profile size =", nrow(x$mp), "\n")
  cat("Window size =", min(x$w), "-", max(x$w), "\n")
  cat("Exclusion zone =", x$ez, "times the windows size\n")

  if (!is.null(x$data)) {
    set <- length(x$data)
    obs <- nrow(x$data[[1]])
    dim <- ncol(x$data[[1]])
    cat(
      "Contains", set, ifelse(set > 1, "sets", "set"), "of data with", obs, "observations and", dim,
      ifelse(dim > 1, "dimensions", "dimension"), "\n"
    )
  }
}

#' Prints a Matrix Profile
#'
#' @param x a TSMP object of class `MatrixProfile`.
#' @param \dots additional arguments ignored.
#' @export
#' @keywords internal
#' @noRd

print.MatrixProfile <- function(x, ...) {
  cat("Matrix Profile\n")
  cat("--------------\n")

  cat("Profile size =", nrow(x$mp), "\n")
  cat("Window size =", x$w, "\n")
  cat("Exclusion zone =", round(x$w * x$ez + vars()$eps), "\n")

  if (!is.null(x$data)) {
    set <- length(x$data)
    obs <- nrow(x$data[[1]])
    dim <- ncol(x$data[[1]])
    cat(
      "Contains", set, ifelse(set > 1, "sets", "set"), "of data with", obs, "observations and", dim,
      ifelse(dim > 1, "dimensions", "dimension"), "\n"
    )
  }
}

#' Prints a Multidimensional Matrix Profile
#'
#' @param x a TSMP object of class `MultiMatrixProfile`.
#' @param \dots additional arguments ignored.
#' @export
#' @keywords internal
#' @noRd


print.MultiMatrixProfile <- function(x, ...) {
  cat("Multidimensional Matrix Profile\n")
  cat("-------------------------------\n")

  cat("Profile size =", nrow(x$mp), "\n")
  cat("Dimensions =", x$n_dim, "\n")
  cat("Window size =", x$w, "\n")
  cat("Exclusion zone =", round(x$w * x$ez + vars()$eps), "\n")
  cat("Must dimensions =", ifelse(is.null(x$must), "None", x$must), "\n")
  cat("Excluded dimensions =", ifelse(is.null(x$exc), "None", x$exc), "\n")

  if (!is.null(x$data)) {
    set <- length(x$data)
    obs <- nrow(x$data[[1]])
    dim <- ncol(x$data[[1]])
    cat(
      "Contains", set, ifelse(set > 1, "sets", "set"), "of data with", obs, "observations and", dim,
      ifelse(dim > 1, "dimensions", "dimension"), "\n"
    )
  }
}

#' Prints a PMP
#'
#' @param x a TSMP object of class `PMP`.
#' @param \dots additional arguments ignored.
#' @export
#' @keywords internal
#' @noRd

print.PMP <- function(x, ...) {
  cat("Pan-Matrix Profile\n")
  cat("------------------\n")

  cat("Number of profiles =", length(x$pmp), "\n")
  cat("Window sizes = from", min(x$w), "to", max(x$w), "\n")
  cat("Exclusion zone =", x$ez, "\n")

  if (!is.null(x$data)) {
    set <- 1
    obs <- length(x$data[[1]])
    cat(
      "Contains", set, ifelse(set > 1, "sets", "set"), "of data with", obs, "observations\n"
    )
  }
}

#' Prints a SiMPle Matrix Profile
#'
#' @param x a TSMP object of class `SimpleMatrixProfile`.
#' @param \dots additional arguments ignored.
#' @export
#' @keywords internal
#' @noRd

print.SimpleMatrixProfile <- function(x, ...) {
  cat("SiMPle Matrix Profile\n")
  cat("---------------------\n")

  cat("Profile size =", nrow(x$mp), "\n")
  cat("Dimensions =", x$n_dim, "\n")
  cat("Window size =", x$w, "\n")
  cat("Exclusion zone =", round(x$w * x$ez + vars()$eps), "\n")

  if (!is.null(x$data)) {
    set <- length(x$data)
    obs <- nrow(x$data[[1]])
    dim <- ncol(x$data[[1]])
    cat(
      "Contains", set, ifelse(set > 1, "sets", "set"), "of data with", obs, "observations and", dim,
      ifelse(dim > 1, "dimensions", "dimension"), "\n"
    )
  }
}

#' Prints a CAC profile
#'
#' @param x a TSMP object of class `ArcCount`.
#' @param \dots additional arguments ignored.
#' @export
#' @keywords internal
#' @noRd
print.ArcCount <- function(x, ...) {
  if (any(class(x) %in% "MatrixProfile")) {
    print.MatrixProfile(x, ...)
  } else if (any(class(x) %in% "MultiMatrixProfile")) {
    print.MultiMatrixProfile(x, ...)
  }

  if (!is.null(x$cac_final)) {
    cat("\nArc Count - Online\n")
    cat("------------------\n")
  } else {
    cat("\nArc Count\n")
    cat("---------\n")
  }

  cat("Profile size =", length(x$cac), "\n")
  cat("Minimum normalized count =", signif(min(x$cac), 2), "at index", which.min(x$cac), "\n")
}

#' Prints a FLOSS
#'
#' @param x a TSMP object of class `Floss`.
#' @param \dots additional arguments ignored.
#' @export
#' @keywords internal
#' @noRd
print.Floss <- function(x, ...) {
  if (any(class(x) %in% "ArcCount")) {
    print.ArcCount(x, ...)
  }

  cat("\nFloss\n")
  cat("-----\n")

  cat("Segments =", length(x$floss), "\n")
  cat("Segmentation indexes =", x$floss, "\n")
  cat("Segmentation thld values =", x$floss_vals, "\n")
}

#' Prints a FLUSS
#'
#' @param x a TSMP object of class `Fluss`.
#' @param \dots additional arguments ignored.
#' @export
#' @keywords internal
#' @noRd
print.Fluss <- function(x, ...) {
  if (any(class(x) %in% "ArcCount")) {
    print.ArcCount(x, ...)
  }

  cat("\nFluss\n")
  cat("-----\n")

  cat("Segments =", length(x$fluss), "\n")
  cat("Segmentation indexes =", x$fluss, "\n")
}

#' Prints a TS Chain
#'
#' @param x a TSMP object of class `Chain`.
#' @param \dots additional arguments ignored.
#' @export
#' @keywords internal
#' @noRd
print.Chain <- function(x, ...) {
  if (any(class(x) %in% "MatrixProfile")) {
    print.MatrixProfile(x, ...)
  } else if (any(class(x) %in% "MultiMatrixProfile")) {
    print.MultiMatrixProfile(x, ...)
  }

  cat("\nChain\n")
  cat("-----\n")

  cat("Chains found =", length(x$chain$chains), "\n")
  cat("Best Chain size =", length(x$chain$best), "\n")
  cat("Best Chain indexes =", x$chain$best, "\n")
}

#' Prints Discords
#'
#' @param x a TSMP object of class `Discord`.
#' @param \dots additional arguments ignored.
#' @export
#' @keywords internal
#' @noRd
print.Discord <- function(x, ...) {
  if (any(class(x) %in% "MatrixProfile")) {
    print.MatrixProfile(x, ...)
  } else if (any(class(x) %in% "MultiMatrixProfile")) {
    print.MultiMatrixProfile(x, ...)
  } else if (any(class(x) %in% "PMP")) {
    print.PMP(x, ...)
  }

  cat("\nDiscord\n")
  cat("-------\n")

  dis_len <- length(x$discord$discord_idx)
  cat("Discords found =", dis_len, "\n")

  indexes <- NULL
  for (i in seq_len(dis_len)) {
    indexes <- paste0(indexes, "[", paste(x$discord$discord_idx[i], collapse = ", "), "] ")
  }

  neighbors <- NULL
  for (i in seq_len(dis_len)) {
    neighbors <- paste0(neighbors, "[", paste(x$discord$discord_neighbor[[i]], collapse = ", "), "] ")
  }

  cat("Discords indexes =", indexes, "\n")
  cat("Discords neighbors =", neighbors, "\n")
}

#' Prints Snippets
#'
#' @param x a TSMP object of class `Snippet`.
#' @param \dots additional arguments ignored.
#' @export
#' @keywords internal
#' @noRd
print.Snippet <- function(x, ...) {
  cat("\nSnippet\n")
  cat("-------\n")

  snip_len <- length(x$snippet_idx)
  cat("Snippets found =", snip_len, "\n")
  cat("Snippets indexes =", x$snippet_idx, "\n")
  cat("Snippets fractions =", sprintf("%1.2f%%", 100 * x$snippet_frac), "\n")
  cat("Snippet size =", x$snippet_size, "\n")
}

#' Prints Motifs
#'
#' @param x a TSMP object of class `Motif`.
#' @param \dots additional arguments ignored.
#' @export
#' @keywords internal
#' @noRd
print.Motif <- function(x, ...) {
  valmod <- FALSE

  if ("Valmod" %in% class(x)) {
    valmod <- TRUE
    print.Valmod(x, ...)
  } else if ("MatrixProfile" %in% class(x)) {
    print.MatrixProfile(x, ...)
  } else if ("MultiMatrixProfile" %in% class(x)) {
    print.MultiMatrixProfile(x, ...)
  } else if (any(class(x) %in% "PMP")) {
    print.PMP(x, ...)
  }

  if (valmod) {
    cat("\nValmod Motif\n")
    cat("------------\n")
  } else {
    cat("\nMotif\n")
    cat("-----\n")
  }

  pairs_len <- length(x$motif$motif_idx)
  cat("Motif pairs found =", pairs_len, "\n")

  indexes <- NULL
  for (i in seq_len(pairs_len)) {
    indexes <- paste0(indexes, "[", paste(x$motif$motif_idx[[i]], collapse = ", "), "] ")
  }

  neighbors <- NULL
  for (i in seq_len(pairs_len)) {
    neighbors <- paste0(neighbors, "[", paste(x$motif$motif_neighbor[[i]], collapse = ", "), "] ")
  }

  cat("Motif pairs indexes =", indexes, "\n")
  cat("Motif pairs neighbors =", neighbors, "\n")

  if (valmod) {
    windows <- NULL
    for (i in seq_len(pairs_len)) {
      windows <- paste0(windows, "[", paste(x$motif$motif_window[[i]], collapse = ", "), "] ")
    }
    cat("Motif pairs windows =", windows, "\n")
  }
}

#' Prints Multidimensional Motifs
#'
#' @param x a TSMP object of class `MultiMotif`.
#' @param \dots additional arguments ignored.
#' @export
#' @keywords internal
#' @noRd
print.MultiMotif <- function(x, ...) {
  if (any(class(x) %in% "MatrixProfile")) {
    print.MatrixProfile(x, ...)
  } else if (any(class(x) %in% "MultiMatrixProfile")) {
    print.MultiMatrixProfile(x, ...)
  }

  cat("\nMultidimensional Motif\n")
  cat("----------------------\n")

  pairs_len <- length(x$motif$motif_idx)
  cat("Motif pairs found =", pairs_len, "\n")

  indexes <- NULL
  for (i in seq_len(pairs_len)) {
    indexes <- paste0(indexes, "[", paste(x$motif$motif_idx[[i]], collapse = ", "), "] ")
  }
  cat("Motif pairs indexes =", indexes, "\n")

  dims <- NULL
  for (i in seq_len(pairs_len)) {
    dims <- paste0(dims, "[", paste(x$motif$motif_dim[[i]], collapse = ", "), "] ")
  }
  cat("Motifs pairs dimensions =", dims, "\n")
}

#' Prints Salient subsequences summary
#'
#' @param x a TSMP object of class `Salient`.
#' @param \dots additional arguments ignored.
#' @export
#' @keywords internal
#' @noRd
print.Salient <- function(x, ...) {
  if (any(class(x) %in% "MatrixProfile")) {
    print.MatrixProfile(x, ...)
  } else if (any(class(x) %in% "MultiMatrixProfile")) {
    print.MultiMatrixProfile(x, ...)
  }

  cat("\nSalient Subsequences\n")
  cat("----------------------\n")

  salient_len <- nrow(x$salient$indexes)

  cat("Subsequences found =", salient_len, "\n")
  cat("Bitsizes tested =", x$salient$bits, "\n")
}
