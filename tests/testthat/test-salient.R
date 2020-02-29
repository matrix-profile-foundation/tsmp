if (!testthat:::on_cran()) {
  context("Testing Salient functions")
  library(tsmp)

  data <- mp_toy_data$data[, 1]
  mp <- tsmp(data, window_size = 30, verbose = 0)
  label_idx <- seq(2, 500, by = 110) # fake data
  res <- NULL


  test_that("Finish", {
    expect_silent(salient_subsequences(mp, verbose = 0))

    if (!testthat:::on_cran()) {
      expect_message(res <<- salient_subsequences(mp, n_bits = c(4, 6, 8), verbose = 2), "Finished")
    } else {
      expect_silent(res <<- salient_subsequences(mp, n_bits = c(4, 6, 8), verbose = 0))
    }
  })

  if (!is.null(res)) {
    test_that("Score", {
      expect_equal(round(sum(res$salient$indexes) / sd(res$salient$indexes), 4), 93.7207)
      expect_equal(round(sum(res$salient$idx_bit_size) / sd(res$salient$idx_bit_size), 2), 195.59)
      expect_equal(tsmp:::get_bitsize(data > 0, 10), 5490)
      expect_equal(sum(tsmp:::discrete_norm(data, 3, max(data), min(data))), 546)
      maxmin <- tsmp:::discrete_norm_pre(as.vector(data), 100)
      expect_equal(round(maxmin$max, 4), 3.3845)
      expect_equal(round(maxmin$min, 4), -3.4308)
      # sort gives slightly different results in windows and linux:
      expect_equal(sum(tsmp:::get_sorted_idx(res$mp, 10) %in% c(36, 408, 37, 407, 35, 200, 9, 199, 10, 406)), 10)
      expect_equal(round(sd(salient_mds(res)), 2), 3.69)
      scr <- salient_score(res, label_idx)
      expect_equal(round(scr$precision, 4), 0.5)
      expect_equal(round(scr$recall, 4), 0.2)
      expect_equal(round(scr$fscore, 4), 0.2857)
    })
  }
}
