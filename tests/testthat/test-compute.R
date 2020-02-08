if (!testthat:::on_cran()) {
  library(tsmp)
  context("Testing if basic functions are ok")

  ts <- mp_toy_data$data[, 1]
  query <- mp_toy_data$data[, 2]
  window <- NULL
  query <- NULL
  sample_pct <- 1.0
  threshold <- 0.98
  n_jobs <- 1L


  test_that("Errors", {
    # windows sizes
    expect_error(compute(ts, -1), "must be >= 4")
    expect_error(compute(ts, query = 1:3), "Must be of length >= 4")
    expect_error(compute(ts, query = c(1:10, NA, 1:10)), "May not contain missing values")
  })
  # single window, no query
  result_single <- compute(ts, windows = 30)


  # single window, and query
  result_query <- compute(ts, windows = 30, query = query)

  # multiple windows
  #
  result_multiple <- compute(ts, windows = c(10:70))
}
