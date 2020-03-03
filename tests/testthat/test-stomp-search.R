if (!testthat:::on_cran()) {
  context("Testing STOMP Search")
  library(tsmp)

  data <- mp_toy_data$data[1:200, 1]
  w <- mp_toy_data$sub_len
  mp <- tsmp(data, window_size = w, mode = "stomp", verbose = 0)
  motifs <- find_motif(mp, n_motifs = 2)

  test_that("Motifs are correct", {
    expect_equal(motifs$motif$motif_idx[[1]], c(33, 102))
    expect_equal(motifs$motif$motif_idx[[2]], c(9, 127))
    expect_equal(length(motifs$motif$motif_neighbor[[1]]), 1)
    expect_equal(motifs$motif$motif_neighbor[[2]], c(148, 77))
    expect_equal(motifs$motif$motif_window[[1]], 30)
  })

  test_that("Errors", {
    mpe <- mp
    mpe$data[[1]] <- table(mp$data[[1]])
    expect_error(find_motif(mpe, n_motifs = 2), "`data` must be `matrix`")
  })

  test_that("Handles other paths", {
    mpe <- mp
    # data as list
    mpe$data[[1]] <- list(mp$data[[1]])
    expect_silent(find_motif(mpe, n_motifs = 2))
    # data as vector
    mpe$data[[1]] <- as.vector(mp$data[[1]])
    expect_silent(find_motif(mpe, n_motifs = 2))
  })
}
