if (!testthat:::on_cran()) {
  context("Testing Time Series Chains")
  library(tsmp)

  w <- 50
  data <- mp_gait_data
  mp <- tsmp(data, window_size = w, exclusion_zone = 1 / 4, verbose = 0)
  res <- find_chains(mp)

  test_that("Find Chains", {
    expect_equal(length(res$chain), 2)
    expect_equal(length(res$chain$chains), 58)
    expect_equal(length(res$chain$best), 6)
    expect_known_hash(res$chain$chains, "d7c3f43152")
  })
}
