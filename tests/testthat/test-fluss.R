if (!testthat:::on_cran()) {
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

  floss <- floss(mp, mp_fluss_data$tilt_abp$data[1001:2000], 1001)

  test_that("Corrected Arc Count", {
    expect_equal(round(mean(cac$cac), 4), 0.9941)
    expect_equal(round(sd(cac$cac), 4), 0.0187)
    expect_equal(round(min(cac$cac), 4), 0.8838)
    expect_equal(max(cac$cac), 1)

    expect_equal(round(mean(floss$cac), 4), 0.8643)
    expect_equal(round(sd(floss$cac), 4), 0.2007)
    expect_equal(round(min(floss$cac), 3), 0)
    expect_equal(max(floss$cac), 1)
    expect_equal(round(mean(floss$cac_final, na.rm = TRUE), 4), 0.9755)
  })

  test_that("Segments found", {
    expect_equal(segments$fluss, c(941, 875, 141))
    expect_equal(floss$floss, 1649)
    expect_equal(round(floss$floss_vals, 3), 0.871)
  })

  test_that("Score", {
    expect_equal(score, 0.259)
  })
}
