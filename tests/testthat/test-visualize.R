if (!testthat:::on_cran()) {
  context("Testing API visualize")
  library(tsmp)
  library(vdiffr)

  ts <- mp_toy_data$data[, 1]

  result_single <- compute(ts, windows = 30)

  plot_test <- function() visualize(result_single)
  plot_test2 <- function() visualize(list(profile = result_single))

  test_that("Plot", {
    expect_doppelganger("visualize", plot_test)
    expect_doppelganger("visualize2", plot_test2)
  })
}
