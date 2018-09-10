if (skip_on_cran()) {
  context("Testing mSTOMP Search")
  library(tsmp)

  data <- mp_toy_data$data[1:200, ]
  w <- mp_toy_data$sub_len
  mp <- tsmp(data, window_size = w, mode = "mstomp", verbose = 0)
  motifs <- find_motif(mp, n_motifs = 2)
}
