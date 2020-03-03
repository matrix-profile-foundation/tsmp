if (!testthat:::on_cran()) {
  context("Testing mSTOMP Search")
  library(tsmp)

  data <- mp_toy_data$data[1:200, ]
  w <- mp_toy_data$sub_len
  mp <- tsmp(data, window_size = w, mode = "mstomp", verbose = 0)
  motifs <- find_motif(mp, n_motifs = 2)
  umotifs <- find_motif(mp, n_motifs = 2, mode = "u")

  test_that("mMotifs are correct", {
    expect_equal(motifs$motif$motif_idx, list(c(45, 108)))
    expect_equal(motifs$motif$motif_dim, list(c(1, 2, 3)))
    expect_null(motifs$motif$motif_neighbor)
    expect_equal(umotifs$motif$motif_idx, list(c(33, 102)))
    expect_equal(umotifs$motif$motif_dim, list(1))
    expect_null(umotifs$motif$motif_neighbor)
  })

  test_that("Errors", {
    mpe <- mp
    mpe$data[[1]] <- table(mp$data[[1]])
    expect_error(find_motif(mpe, n_motifs = 2), "`data` must be `matrix`")
  })

  test_that("Handles other paths", {
    mpe <- mp
    # data as list
    mpe$data[[1]] <- list(mp$data[[1]][, 1], mp$data[[1]][, 2], mp$data[[1]][, 3])
    expect_silent(find_motif(mpe, n_motifs = 2))
    # data as vector
    mpe$data[[1]] <- as.vector(mp$data[[1]])
    expect_warning(find_motif(mpe, n_motifs = 2), "`data` dimensions are different")
  })
}
