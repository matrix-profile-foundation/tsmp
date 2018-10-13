#' Prints a Matrix Profile
#'
#' @param x a TSMP object of class `MatrixProfile`.
#' @param ... additional arguments ignored.
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
#' @param ... additional arguments ignored.
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

#' Prints a CAC profile
#'
#' @param x a TSMP object of class `ArcCount`.
#' @param ... additional arguments ignored.
#' @export
#' @keywords internal
#' @noRd
print.ArcCount <- function(x, ...) {
  if (any(class(x) %in% "MatrixProfile")) {
    print.MatrixProfile(x, ...)
  } else if (any(class(x) %in% "MultiMatrixProfile")) {
    print.MultiMatrixProfile(x, ...)
  }

  cat("\nArc Count\n")
  cat("---------\n")

  cat("Profile size =", length(x$cac), "\n")
  cat("Minimum normalized count =", signif(min(x$cac), 2), "at index", which.min(x$cac), "\n")
}

#' Prints a FLUSS
#'
#' @param x a TSMP object of class `Fluss`.
#' @param ... additional arguments ignored.
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
#' @param ... additional arguments ignored.
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

  cat("Chains founded =", length(x$chain$chains), "\n")
  cat("Best Chain size =", length(x$chain$best), "\n")
  cat("Best Chain indexes =", x$chain$best, "\n")
}

#' Prints Motifs
#'
#' @param x a TSMP object of class `Motif`.
#' @param ... additional arguments ignored.
#' @export
#' @keywords internal
#' @noRd
print.Motif <- function(x, ...) {
  if (any(class(x) %in% "MatrixProfile")) {
    print.MatrixProfile(x, ...)
  } else if (any(class(x) %in% "MultiMatrixProfile")) {
    print.MultiMatrixProfile(x, ...)
  }

  cat("\nMotif\n")
  cat("-----\n")

  pairs_len <- length(x$motif$motif_idx)
  cat("Motif pairs founded =", pairs_len, "\n")

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
}

#' Prints Multidimensional Motifs
#'
#' @param x a TSMP object of class `MultiMotif`.
#' @param ... additional arguments ignored.
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

  motifs_len <- length(x$motif$motif_idx)

  dims <- NULL
  for (i in seq_len(motifs_len)) {
    dims <- paste0(dims, "[", paste(x$motif$motif_dim[[i]], collapse = ", "), "] ")
  }
  cat("Motifs founded =", motifs_len, "\n")
  cat("Motifs indexes =", x$motif$motif_idx, "\n")
  cat("Motifs dimensions =", dims, "\n")
}

#' Prints Salient subsequences summary
#'
#' @param x a TSMP object of class `Salient`.
#' @param ... additional arguments ignored.
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

  cat("Subsequences founded =", salient_len, "\n")
  cat("Bitsizes tested =", x$salient$bits, "\n")
}
