#' Title
#' @export
plot.MatrixProfile <- function(.mp, ...) {
  message("DEBUG: calling ", match.call()[[1]])
  def_par <- graphics::par(no.readonly = TRUE)

  graphics::par(def_par)
}
#' Title
#' @export
plot.MultiMatrixProfile <- function(.mp, ...) {
  message("DEBUG: calling ", match.call()[[1]])
  def_par <- graphics::par(no.readonly = TRUE)

  graphics::par(def_par)
}
#' Title
#' @export
plot.Arcs <- function(.mp, pairs, alpha = NULL, quality = 30, lwd = 15, col = c("blue", "orange"),
                      main = "Arc Plot", ylab = "", xlab = "Profile Index", ...) {
  message("DEBUG: calling ", match.call()[[1]])
  def_par <- graphics::par(no.readonly = TRUE)

  xmin <- min(pairs)
  xmax <- max(pairs)
  plot_size <- (xmax - xmin)

  z_seq <- seq(0, base::pi, length.out = quality)
  xlim <- c(xmin - lwd, xmax + lwd)
  ylim <- c(0, plot_size / 2)

  if (is.null(alpha)) {
    alpha <- min(0.5, max(10 / nrow(pairs), 0.03))
  }

  arccolr <- adjustcolor(col, alpha.f = alpha)
  if (length(col) > 1) {
    arccoll <- adjustcolor(col[2], alpha.f = alpha)
  } else {
    arccoll <- adjustcolor(col, alpha.f = alpha)
  }

  # blank plot
  graphics::plot(0.5, 0.5,
    type = "n", main = main, xlab = xlab, ylab = ylab,
    xlim = xlim, ylim = ylim, yaxt = "n", ...
  )

  for (i in seq_len(nrow(pairs))) {
    if (pairs[i, 1] > pairs[i, 2]) {
      arccol <- arccoll
    } else {
      arccol <- arccolr
    }

    x1 <- min(pairs[i, 1], pairs[i, 2])
    x2 <- max(pairs[i, 1], pairs[i, 2])
    center <- (x1 - x2) / 2 + x2
    radius <- (x2 - x1) / 2
    x_seq <- center + radius * cos(z_seq)
    y_seq <- radius * sin(z_seq)
    graphics::lines(x_seq, y_seq,
      col = arccol, lwd = lwd, lty = 1, lend = 1
    )
  }

  graphics::legend(xmin, plot_size / 2,
    legend = c("Right", "Left"),
    col = adjustcolor(col, alpha.f = 0.5), lty = 1, cex = 0.8, lwd = 5
  )

  graphics::par(def_par)
}
#' Title
#' @export
plot.ArcCount <- function(.mp, ...) {
  message("DEBUG: calling ", match.call()[[1]])
  def_par <- graphics::par(no.readonly = TRUE)

  graphics::par(def_par)
}
#' Title
#' @export
plot.Chain <- function(.mp, ...) {
  message("DEBUG: calling ", match.call()[[1]])
  def_par <- graphics::par(no.readonly = TRUE)
  chain_size <- length(.mp$chain$best)
  pairs <- matrix(Inf, round(chain_size / 2) * 2, 2)

  for (i in seq_len(chain_size)) {
    if (i == chain_size) {
      break
    } else {
      pairs[i, 1] <- .mp$chain$best[i]
      pairs[i, 2] <- .mp$chain$best[i + 1]
    }
  }

  plot.Arcs(.mp, pairs)

  # motifs <- .mp$motif$motif_idx
  # n_motifs <- length(.mp$motif$motif_idx)
  # neighbors <- .mp$motif$motif_neighbor
  # matrix_profile_size <- nrow(.mp$mp)
  #
  # # layout: matrix profile on top, motifs below.
  # layout(matrix(
  #   c(rep(1, cols), (seq_len(ceiling(n_motifs / cols) * cols) + 1)),
  #   ceiling(n_motifs / cols) + 1,
  #   cols,
  #   byrow = TRUE
  # ))
  # # plot matrix profile
  # graphics::par(oma = c(1, 1, 4, 0), cex.lab = 1.5)
  # graphics::plot(.mp$mp, type = "l", main = "Matrix Profile", xlab = "index", ylab = "distance")
  # graphics::mtext("MOTIF Discover", line = 4, font = 2, cex = 1.5)
  # abline(v = unlist(motifs), col = rep(1:n_motifs, each = 2), lwd = 2)
  # # plot motifs
  # for (i in 1:n_motifs) {
  #   motif1 <- znorm(data[motifs[[i]][1]:min((motifs[[i]][1] + .mp$w - 1), matrix_profile_size)])
  #   motif2 <- znorm(data[motifs[[i]][2]:min((motifs[[i]][2] + .mp$w - 1), matrix_profile_size)])
  #
  #   # blank plot
  #   graphics::plot(0.5, 0.5,
  #                  type = "n", main = paste("Motif", i), xlab = "length", ylab = "normalized data",
  #                  xlim = c(0, length(motif1)), ylim = c(min(motif1), max(motif1))
  #   )
  #
  #   for (j in seq_len(length(neighbors[[i]]))) {
  #     neigh <- znorm(data[neighbors[[i]][j]:min((neighbors[[i]][j] + .mp$w - 1), matrix_profile_size)])
  #     graphics::lines(neigh, col = "gray70")
  #   }
  #
  #   graphics::lines(motif2, col = "black")
  #   graphics::lines(motif1, col = i, lwd = 2)
  # }

  graphics::par(def_par)
}
#' Title
#' @export
plot.Motif <- function(.mp, data, cols = 3, ...) {
  message("DEBUG: calling ", match.call()[[1]])
  def_par <- graphics::par(no.readonly = TRUE)

  if (missing(data) && !is.null(.mp$data)) {
    data <- .mp$data[[1]]
  } else {
    is.null(data) # check data presence before plotting anything
  }

  motifs <- .mp$motif$motif_idx
  n_motifs <- length(.mp$motif$motif_idx)
  neighbors <- .mp$motif$motif_neighbor
  matrix_profile_size <- nrow(.mp$mp)

  # layout: matrix profile on top, motifs below.
  layout(matrix(
    c(rep(1, cols), (seq_len(ceiling(n_motifs / cols) * cols) + 1)),
    ceiling(n_motifs / cols) + 1,
    cols,
    byrow = TRUE
  ))
  # plot matrix profile
  graphics::par(oma = c(1, 1, 4, 0), cex.lab = 1.5)
  graphics::plot(.mp$mp, type = "l", main = "Matrix Profile", xlab = "index", ylab = "distance")
  graphics::mtext("MOTIF Discover", line = 4, font = 2, cex = 1.5)
  abline(v = unlist(motifs), col = rep(1:n_motifs, each = 2), lwd = 2)
  # plot motifs
  for (i in 1:n_motifs) {
    motif1 <- znorm(data[motifs[[i]][1]:min((motifs[[i]][1] + .mp$w - 1), matrix_profile_size)])
    motif2 <- znorm(data[motifs[[i]][2]:min((motifs[[i]][2] + .mp$w - 1), matrix_profile_size)])

    # blank plot
    graphics::plot(0.5, 0.5,
      type = "n", main = paste("Motif", i), xlab = "length", ylab = "normalized data",
      xlim = c(0, length(motif1)), ylim = c(min(motif1), max(motif1))
    )

    for (j in seq_len(length(neighbors[[i]]))) {
      neigh <- znorm(data[neighbors[[i]][j]:min((neighbors[[i]][j] + .mp$w - 1), matrix_profile_size)])
      graphics::lines(neigh, col = "gray70")
    }

    graphics::lines(motif2, col = "black")
    graphics::lines(motif1, col = i, lwd = 2)
  }

  graphics::par(def_par)
}
#' Title
#' @export
plot.MultiMotif <- function(.mp, data, ...) {
  message("DEBUG: calling ", match.call()[[1]])
  def_par <- graphics::par(no.readonly = TRUE)

  if (missing(data) && !is.null(.mp$data)) {
    data <- .mp$data[[1]]
  } else {
    is.null(data) # check data presence before plotting anything
  }

  graphics::par(def_par)
}
