if (!testthat:::on_cran()) {
  context("Testing API visualize")
  library(tsmp)
  library(vdiffr)

  ts <- mp_toy_data$data[, 1]

  result_single <- compute(ts, windows = 30)

  plot_test <- function() visualize(result_single)

  test_that("Plot", {
    expect_doppelganger("plot", plot_test)
  })
}
