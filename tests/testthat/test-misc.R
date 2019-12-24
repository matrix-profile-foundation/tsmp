if (!testthat:::on_cran()) {
  context("Testing Class Functions")
  library(tsmp)

  data <- mp_toy_data$data[1:100, ]
  w <- 10
  nseg <- 3
  mp <- tsmp(data[, 1], window_size = w, verbose = 0)
  mmp <- tsmp(data, mode = "mstomp", window_size = w, verbose = 0)
  cac <- fluss_cac(mp)
  segments <- fluss_extract(cac, nseg)
  motifs <- find_motif(mp)
  mmotifs <- find_motif(mmp)
  chains <- find_chains(mp)

  test_that("as. functions", {
    expect_equal(class(as.matrixprofile(segments))[1], "MatrixProfile")
    expect_equal(class(as.arccount(segments))[1], "ArcCount")
    expect_equal(class(as.fluss(as.arccount(segments)))[1], "Fluss")
    expect_equal(class(as.chain(as.matrixprofile(chains)))[1], "Chain")
    expect_equal(class(as.motif(as.matrixprofile(motifs)))[1], "Motif")
    expect_equal(class(as.multimatrixprofile(mmotifs))[1], "MultiMatrixProfile")
    expect_equal(class(as.multimotif(as.multimatrixprofile(mmotifs)))[1], "MultiMotif")

    expect_error(as.matrixprofile(mmotifs), "cannot be")
    expect_error(as.arccount(mmotifs), "cannot be")
    expect_error(as.fluss(mmotifs), "cannot be")
    expect_error(as.chain(mmotifs), "cannot be")
    expect_error(as.motif(mmotifs), "cannot be")
    expect_error(as.multimatrixprofile(segments), "cannot be")
    expect_error(as.multimotif(segments), "cannot be")
    expect_error(as.salient(segments), "cannot be")
  })
}
