if (skip_on_cran()) {
  context("Testing Subset functions")
  library(tsmp)

  cac <- fluss_cac(test_mp)
  segments <- fluss_extract(cac, nseg)

  test_that("Corrected Arc Count", {
    expect_equal(round(mean(cac$cac), 4), 0.5296)
    expect_equal(round(sd(cac$cac), 3), 0.404)
    expect_equal(round(min(cac$cac), 3), 0)
    expect_equal(max(cac$cac), 1)
  })

  test_that("Head Arc Count", {
    h_cac <- head(cac, 20000)
    expect_equal(round(mean(h_cac$cac), 4), 0.6501)
    expect_equal(round(sd(h_cac$cac), 3), 0.328)
    expect_equal(round(min(h_cac$cac), 3), 0.152)
    expect_equal(max(h_cac$cac), 1)
  })

  test_that("Tail Arc Count", {
    t_cac <- tail(cac, 20000)
    expect_equal(round(mean(t_cac$cac), 4), 0.5327)
    expect_equal(round(sd(t_cac$cac), 3), 0.395)
    expect_equal(round(min(t_cac$cac), 3), 0)
    expect_equal(max(t_cac$cac), 1)
  })

  test_that("Segments found", {
    expect_equal(segments$fluss, c(24907, 26217, 27643))
  })

  test_that("Head Segments found", {
    h_segment <- head(segments, 20000)
    expect_equal(h_segment$fluss, c(7150, 8445, 5852))
  })

  test_that("Tail Segments found", {
    t_segment <- tail(segments, 20000)
    expect_equal(t_segment$fluss, c(4664, 6217, 7643))
  })
}
