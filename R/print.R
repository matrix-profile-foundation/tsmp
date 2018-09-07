#' Prints a Matrix Profile
#' @export
#' @keywords internal
#' @noRd
print.MatrixProfile <- function(.mp, ...) {
  message("DEBUG: calling ", match.call()[[1]])

  cat("Matrix Profile\n")
  cat("--------------\n")

  cat("Profile size =", nrow(.mp$mp), "\n")
  cat("Window size =", .mp$w, "\n")
  cat("Exclusion zone =", round(.mp$w * .mp$ez + vars()$eps), "\n")

  if (!is.null(.mp$data)) {
    set <- length(.mp$data)
    obs <- nrow(.mp$data[[1]])
    dim <- ncol(.mp$data[[1]])
    cat(
      "Contains", set, ifelse(set > 1, "sets", "set"), "of data with", obs, "observations and", dim,
      ifelse(dim > 1, "dimensions", "dimension"), "\n"
    )
  }
}

#' Prints a Multidimensional Matrix Profile
#' @export
#' @keywords internal
#' @noRd
print.MultiMatrixProfile <- function(.mp, ...) {
  message("DEBUG: calling ", match.call()[[1]])

  cat("Multidimensional Matrix Profile\n")
  cat("-------------------------------\n")

  cat("Profile size =", nrow(.mp$mp), "\n")
  cat("Dimensions =", .mp$n_dim, "\n")
  cat("Window size =", .mp$w, "\n")
  cat("Exclusion zone =", round(.mp$w * .mp$ez + vars()$eps), "\n")
  cat("Must dimensions =", ifelse(is.null(.mp$must), "None", .mp$must), "\n")
  cat("Excluded dimensions =", ifelse(is.null(.mp$exc), "None", .mp$exc), "\n")

  if (!is.null(.mp$data)) {
    set <- length(.mp$data)
    obs <- nrow(.mp$data[[1]])
    dim <- ncol(.mp$data[[1]])
    cat(
      "Contains", set, ifelse(set > 1, "sets", "set"), "of data with", obs, "observations and", dim,
      ifelse(dim > 1, "dimensions", "dimension"), "\n"
    )
  }
}

#' Prints a CAC profile
#' @export
#' @keywords internal
#' @noRd
print.ArcCount <- function(.mp, ...) {
  message("DEBUG: calling ", match.call()[[1]])

  if (any(class(.mp) %in% c("MatrixProfile"))) {
    print.MatrixProfile(.mp, ...)
  } else if (any(class(.mp) %in% c("MultiMatrixProfile"))) {
    print.MultiMatrixProfile(.mp, ...)
  }

  cat("\nArc Count\n")
  cat("---------\n")

  cat("Profile size =", length(.mp$cac), "\n")
  cat("Minimum normalized count =", signif(min(.mp$cac), 2), "at index", which.min(.mp$cac), "\n")
}

#' Prints a FLUSS
#' @export
#' @keywords internal
#' @noRd
print.Fluss <- function(.mp, ...) {
  message("DEBUG: calling ", match.call()[[1]])

  if (any(class(.mp) %in% c("ArcCount"))) {
    print.ArcCount(.mp, ...)
  }

  cat("\nFluss\n")
  cat("-----\n")

  cat("Segments =", length(.mp$fluss), "\n")
  cat("Segmentation indexes =", .mp$fluss, "\n")
}

#' Prints a TS Chain
#' @export
#' @keywords internal
#' @noRd
print.Chain <- function(.mp, ...) {
  message("DEBUG: calling ", match.call()[[1]])

  if (any(class(.mp) %in% c("MatrixProfile"))) {
    print.MatrixProfile(.mp, ...)
  } else if (any(class(.mp) %in% c("MultiMatrixProfile"))) {
    print.MultiMatrixProfile(.mp, ...)
  }

  cat("\nChain\n")
  cat("-----\n")

  cat("Chains founded =", length(.mp$chain$chains), "\n")
  cat("Best Chain size =", length(.mp$chain$best), "\n")
  cat("Best Chain indexes =", .mp$chain$best, "\n")
}

#' Prints Motifs
#' @export
#' @keywords internal
#' @noRd
print.Motif <- function(.mp, ...) {
  message("DEBUG: calling ", match.call()[[1]])

  if (any(class(.mp) %in% c("MatrixProfile"))) {
    print.MatrixProfile(.mp, ...)
  } else if (any(class(.mp) %in% c("MultiMatrixProfile"))) {
    print.MultiMatrixProfile(.mp, ...)
  }

  cat("\nMotif\n")
  cat("-----\n")

  pairs_len <- length(.mp$motif$motif_idx)
  cat("Motif pairs founded =", pairs_len, "\n")

  indexes <- NULL
  for (i in seq_len(pairs_len)) {
    indexes <- paste0(indexes, "[", paste(.mp$motif$motif_idx[[i]], collapse = ", "), "] ")
  }

  neighbors <- NULL
  for (i in seq_len(pairs_len)) {
    neighbors <- paste0(neighbors, "[", paste(.mp$motif$motif_neighbor[[i]], collapse = ", "), "] ")
  }

  cat("Motif pairs indexes =", indexes, "\n")
  cat("Motif pairs neighbors =", neighbors, "\n")
}

#' Prints Multidimensional Motifs
#' @export
#' @keywords internal
#' @noRd
print.MultiMotif <- function(.mp, ...) {
  message("DEBUG: calling ", match.call()[[1]])

  if (any(class(.mp) %in% c("MatrixProfile"))) {
    print.MatrixProfile(.mp, ...)
  } else if (any(class(.mp) %in% c("MultiMatrixProfile"))) {
    print.MultiMatrixProfile(.mp, ...)
  }

  cat("\nMultidimensional Motif\n")
  cat("----------------------\n")

  motifs_len <- length(.mp$motif$motif_idx)

  dims <- NULL
  for (i in seq_len(motifs_len)) {
    dims <- paste0(dims, "[", paste(.mp$motif$motif_dim[[i]], collapse = ", "), "] ")
  }
  cat("Motifs founded =", motifs_len, "\n")
  cat("Motifs indexes =", .mp$motif$motif_idx, "\n")
  cat("Motifs dimensions =", dims, "\n")
}
