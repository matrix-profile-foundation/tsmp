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
  discords <- find_discord(mp)
  mmotifs <- find_motif(mmp)
  chains <- find_chains(mp)
  corr <- ed_corr(mp$mp, w)
  norm <- normalize(data[, 1])
  paa_t <- paa(data[, 1], 2)
  paa_i <- ipaa(paa_t, 2)

  test_that("as. functions", {
    expect_equal(class(as.matrixprofile(segments))[1], "MatrixProfile")
    expect_equal(class(as.arccount(segments))[1], "ArcCount")
    expect_equal(class(as.fluss(as.arccount(segments)))[1], "Fluss")
    expect_equal(class(as.chain(as.matrixprofile(chains)))[1], "Chain")
    expect_equal(class(as.discord(as.matrixprofile(discords)))[1], "Discord")
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

  test_that("misc", {
    expect_silent(set_data(mp, data[, 1]))
    expect_silent(remove_class(motifs, "Motif"))
    expect_equal(as.vector(data[, 1]), as.vector(get_data(mp)))
    expect_equal(round(sum(corr) / sd(corr), 3), 695.805)
    expect_equal(round(sum(norm) / sd(norm), 3), 168.874)
    expect_equal(round(sum(paa_t) / sd(paa_t), 3), 155.269)
    expect_equal(round(sum(paa_i) / sd(paa_i), 3), 312.118)
  })
}
