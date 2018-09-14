if (skip_on_cran()) {
  context("Testing if basic functions are ok")
  library(tsmp)
  w <- 30
  ref_data <- mp_toy_data$data[, 1]
  query_data <- mp_toy_data$data[, 1]
  d_size <- length(ref_data)
  q_size <- length(query_data)

  test_that("Plot", {
    expect_p
  })

  test_that("Errors", {
    # big window size
    expect_error(fast_movsd(mp_toy_data$data[, 1], 1), regexp = "must be at least 2")
    expect_error(fast_movsd(mp_toy_data$data[1:100, 1], 500), regexp = "is too large")
    expect_error(expect_message(beep(audio::close.audioInstance(99)), "Failed"))
    expect_error(tsmp:::diff2(data.frame(1:10), as.matrix(10:1)), regexp = "matrices")
    expect_error(tsmp:::diff2(as.matrix(1:10), matrix(10:1, ncol = 2)), regexp = "columns")
  })
}
