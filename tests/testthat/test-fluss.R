if (skip_on_cran()) {
  context("Testing FLUSS functions")
  library(tsmp)

  data <- mp_fluss_data$tilt_abp$data[1:1000]
  w <- 10
  truth <- 400
  nseg <- 3
  mp <- tsmp(data, window_size = w, verbose = 0)
  cac <- fluss_cac(mp)
  segments <- fluss_extract(cac, nseg)
  score <- fluss_score(truth, segments$fluss, length(data))

  test_that("Corrected Arc Count", {
    expect_equal(round(mean(cac$cac), 4), 0.9941)
    expect_equal(round(sd(cac$cac), 4), 0.0187)
    expect_equal(round(min(cac$cac), 4), 0.8838)
    expect_equal(max(cac$cac), 1)
  })

  test_that("Segments found", {
    expect_equal(segments$fluss, c(941, 875, 141))
  })

  test_that("Score", {
    expect_equal(score, 0.259)
  })
}
