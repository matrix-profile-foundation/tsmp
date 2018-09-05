plot_motif <- function(motif_idxs) {
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
