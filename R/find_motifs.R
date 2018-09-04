#' Find MOTIFs and Plot them
#'
#' @param matrix_profile
#' @param profile_index
#' @param window_size
#' @param data
#' @param exclusion_zone
#' @param radius
#' @param n_motifs
#' @param plot
#' @param cols
#'
#' @return
#' @export
#'
#' @examples
find_motifs <- function(matrix_profile, profile_index,
                        window_size, data,
                        exclusion_zone = 1/2, radius = 3, n_motifs = 3, plot = TRUE, cols = 3) {

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

  exclusion_zone <- floor(window_size * exclusion_zone)
  data_size <- nrow(data)
  matrix_profile_size <- length(matrix_profile)
  motif_idxs <- list(motifs = list(NULL), neighbors = list(NULL))
  mp <- matrix_profile

  nn_pre <- mass_pre(data, data_size, window_size = window_size)

  for (i in 1:n_motifs) {
    min_idx <- which.min(matrix_profile)
    motif_distance <- matrix_profile[min_idx]
    motif_distance <- motif_distance^2
    motif_idxs[[1]][[i]] <- sort(c(min_idx, profile_index[min_idx]))
    motif_idx <- motif_idxs[[1]][[i]][1]

    query <- data[motif_idx:(motif_idx + window_size - 1)]

    distance_profile <- mass(
      nn_pre$data_fft, query, data_size, window_size, nn_pre$data_mean, nn_pre$data_sd,
      nn_pre$data_mean[motif_idx], nn_pre$data_sd[motif_idx]
    )

    distance_profile <- Re(distance_profile$distance_profile)
    distance_profile[distance_profile > motif_distance * radius] <- Inf
    motif_zone_start <- max(1, motif_idx - exclusion_zone)
    motif_zone_end <- min(matrix_profile_size, motif_idx + exclusion_zone)
    distance_profile[motif_zone_start:motif_zone_end] <- Inf
    motif_idx <- motif_idxs[[1]][[i]][2]
    motif_zone_start <- max(1, motif_idx - exclusion_zone)
    motif_zone_end <- min(matrix_profile_size, motif_idx + exclusion_zone)
    distance_profile[motif_zone_start:motif_zone_end] <- Inf
    st <- sort(distance_profile, index.return = TRUE)
    distance_order <- st$x
    distance_idx_order <- st$ix

    motif_neighbor <- vector(mode = "numeric")

    for (j in 1:10) {
      if (is.infinite(distance_order[1]) || length(distance_order) < j) {
        break
      }
      motif_neighbor[j] <- distance_idx_order[1]
      distance_order <- distance_order[2:length(distance_order)]
      distance_idx_order <- distance_idx_order[2:length(distance_idx_order)]
      distance_order <- distance_order[!(abs(distance_idx_order - motif_neighbor[j]) < exclusion_zone)]
      distance_idx_order <- distance_idx_order[!(abs(distance_idx_order - motif_neighbor[j]) < exclusion_zone)]
    }

    motif_neighbor <- motif_neighbor[motif_neighbor != 0]
    motif_idxs[[2]][[i]] <- motif_neighbor

    remove_idx <- c(motif_idxs[[1]][[i]], motif_idxs[[2]][[i]])

    for (j in 1:length(remove_idx)) {
      remove_zone_start <- max(1, remove_idx[j] - exclusion_zone)
      remove_zone_end <- min(matrix_profile_size, remove_idx[j] + exclusion_zone)
      matrix_profile[remove_zone_start:remove_zone_end] <- Inf
    }
  }

  if (plot == TRUE) {
    def_par <- par(no.readonly = TRUE)
    # layout: matrix profile on top, motifs below.
    layout(matrix(c(rep(1, cols), (seq_len(ceiling(n_motifs / cols) * cols) + 1)),
      ceiling(n_motifs / cols) + 1,
      cols,
      byrow = TRUE
    ))
    # plot matrix profile
    par(oma = c(1, 1, 4, 0), cex.lab = 1.5)
    plot(mp, type = "l", main = "Matrix Profile", xlab = "index", ylab = "distance")
    mtext("MOTIF Discover", line = 4, font = 2, cex = 1.5)
    abline(v = unlist(motif_idxs[[1]]), col = rep(1:n_motifs, each = 2), lwd = 2)
    # plot motifs
    for (i in 1:n_motifs) {
      motif1 <- znorm(data[motif_idxs[[1]][[i]][1]:min((motif_idxs[[1]][[i]][1] +
        window_size - 1), matrix_profile_size)])
      motif2 <- znorm(data[motif_idxs[[1]][[i]][2]:min((motif_idxs[[1]][[i]][2] +
        window_size - 1), matrix_profile_size)])

      # blank plot
      plot(0.5, 0.5,
        type = "n", main = paste("Motif", i), xlab = "length", ylab = "normalized data",
        xlim = c(0, length(motif1)), ylim = c(min(motif1), max(motif1))
      )

      for (j in seq_len(length(motif_idxs[[2]][[i]]))) {
        neigh <- znorm(data[motif_idxs[[2]][[i]][j]:min((motif_idxs[[2]][[i]][j] +
          window_size - 1), matrix_profile_size)])
        lines(neigh, col = "gray70")
      }

      lines(motif2, col = "black")
      lines(motif1, col = i, lwd = 2)
    }
    par(def_par)
  }

  return(motif_idxs)
}
