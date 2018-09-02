#' Find MOTIFs and Plot them
#'
#' @param matrix.profile
#' @param profile.index
#' @param window.size
#' @param data
#' @param exclusion.zone
#' @param radius
#' @param n.motifs
#' @param plot
#' @param cols
#'
#' @return
#' @export
#'
#' @examples
find.motifs <- function(matrix.profile, profile.index,
                        window.size, data,
                        exclusion.zone = 1/2, radius = 3, n.motifs = 3, plot = TRUE, cols = 3) {

  # transform data into matrix
  if (is.vector(data)) {
    data <- as.matrix(data)
  }
  else if (is.matrix(data)) {
    if (ncol(data) > nrow(data)) {
      data <- t(data)
    }
  } else {
    stop("Error: Unknown type of data. Must be: a column matrix or a vector.", call. = FALSE)
  }

  exclusion.zone <- floor(window.size * exclusion.zone)
  data.size <- nrow(data)
  matrix.profile.size <- length(matrix.profile)
  motif.idxs <- list(motifs = list(NULL), neighbors = list(NULL))
  mp <- matrix.profile

  nn.pre <- mass.pre(data, data.size, window.size = window.size)

  for (i in 1:n.motifs) {
    min.idx <- which.min(matrix.profile)
    motif.distance <- matrix.profile[min.idx]
    motif.distance <- motif.distance^2
    motif.idxs[[1]][[i]] <- sort(c(min.idx, profile.index[min.idx]))
    motif.idx <- motif.idxs[[1]][[i]][1]

    query <- data[motif.idx:(motif.idx + window.size - 1)]

    distance.profile <- mass(
      nn.pre$data.fft, query, data.size, window.size, nn.pre$data.mean, nn.pre$data.sd,
      nn.pre$data.mean[motif.idx], nn.pre$data.sd[motif.idx]
    )

    distance.profile <- Re(distance.profile$distance.profile)
    distance.profile[distance.profile > motif.distance * radius] <- Inf
    motif.zone.start <- max(1, motif.idx - exclusion.zone)
    motif.zone.end <- min(matrix.profile.size, motif.idx + exclusion.zone)
    distance.profile[motif.zone.start:motif.zone.end] <- Inf
    motif.idx <- motif.idxs[[1]][[i]][2]
    motif.zone.start <- max(1, motif.idx - exclusion.zone)
    motif.zone.end <- min(matrix.profile.size, motif.idx + exclusion.zone)
    distance.profile[motif.zone.start:motif.zone.end] <- Inf
    st <- sort(distance.profile, index.return = TRUE)
    distance.order <- st$x
    distance.idx.order <- st$ix

    motif.neighbor <- vector(mode = "numeric")

    for (j in 1:10) {
      if (is.infinite(distance.order[1]) || length(distance.order) < j) {
        break
      }
      motif.neighbor[j] <- distance.idx.order[1]
      distance.order <- distance.order[2:length(distance.order)]
      distance.idx.order <- distance.idx.order[2:length(distance.idx.order)]
      distance.order <- distance.order[!(abs(distance.idx.order - motif.neighbor[j]) < exclusion.zone)]
      distance.idx.order <- distance.idx.order[!(abs(distance.idx.order - motif.neighbor[j]) < exclusion.zone)]
    }

    motif.neighbor <- motif.neighbor[motif.neighbor != 0]
    motif.idxs[[2]][[i]] <- motif.neighbor

    remove.idx <- c(motif.idxs[[1]][[i]], motif.idxs[[2]][[i]])

    for (j in 1:length(remove.idx)) {
      remove.zone.start <- max(1, remove.idx[j] - exclusion.zone)
      remove.zone.end <- min(matrix.profile.size, remove.idx[j] + exclusion.zone)
      matrix.profile[remove.zone.start:remove.zone.end] <- Inf
    }
  }

  if (plot == TRUE) {
    def.par <- par(no.readonly = TRUE)
    # layout: matrix profile on top, motifs below.
    layout(matrix(c(rep(1, cols), (seq_len(ceiling(n.motifs / cols) * cols) + 1)),
      ceiling(n.motifs / cols) + 1,
      cols,
      byrow = TRUE
    ))
    # plot matrix profile
    par(oma = c(1, 1, 4, 0), cex.lab = 1.5)
    plot(mp, type = "l", main = "Matrix Profile", xlab = "index", ylab = "distance")
    mtext("MOTIF Discover", line = 4, font = 2, cex = 1.5)
    abline(v = unlist(motif.idxs[[1]]), col = rep(1:n.motifs, each = 2), lwd = 2)
    # plot motifs
    for (i in 1:n.motifs) {
      motif1 <- znorm(data[motif.idxs[[1]][[i]][1]:min((motif.idxs[[1]][[i]][1] +
        window.size - 1), matrix.profile.size)])
      motif2 <- znorm(data[motif.idxs[[1]][[i]][2]:min((motif.idxs[[1]][[i]][2] +
        window.size - 1), matrix.profile.size)])

      # blank plot
      plot(0.5, 0.5,
        type = "n", main = paste("Motif", i), xlab = "length", ylab = "normalized data",
        xlim = c(0, length(motif1)), ylim = c(min(motif1), max(motif1))
      )

      for (j in seq_len(length(motif.idxs[[2]][[i]]))) {
        neigh <- znorm(data[motif.idxs[[2]][[i]][j]:min((motif.idxs[[2]][[i]][j] +
          window.size - 1), matrix.profile.size)])
        lines(neigh, col = "gray70")
      }

      lines(motif2, col = "black")
      lines(motif1, col = i, lwd = 2)
    }
    par(def.par)
  }

  return(motif.idxs)
}
