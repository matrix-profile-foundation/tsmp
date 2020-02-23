if (!testthat:::on_cran()) {
  context("Testing API analyze")
  library(tsmp)

  ts <- mp_toy_data$data[, 1]
  query <- mp_toy_data$data[, 2]
  window <- NULL
  sample_pct <- 1.0
  threshold <- 0.98
  n_jobs <- 1L

  test_that("Errors", {
    # windows sizes
    expect_error(analyze(ts, -1), "must be >= 4")
    expect_error(analyze(ts, query = 1:3), "Must be of length >= 4")
    expect_error(analyze(ts, query = c(1:10, NA, 1:10)), "May not contain missing values")
  })

  result_single <- analyze(ts, windows = 30)

  test_that("Result Single", {
    expect_named(result_single)
    expect_type(result_single, "list")
    expect_s3_class(result_single, "MatrixProfile")
    expect_length(result_single, 12)
    expect_false(attr(result_single, "join"))
    expect_equal(round(mean(result_single$mp), 3), 2.817)
    expect_equal(round(sd(result_single$mp), 4), 0.8975)
    expect_equal(round(mean(result_single$pi), 4), 241.0192)
    expect_equal(round(sd(result_single$pi), 4), 157.7044)
  })


  # single window, and query
  result_query <- analyze(ts, windows = 30, query = query)

  test_that("Result Query", {
    expect_named(result_query)
    expect_type(result_query, "list")
    expect_s3_class(result_query, "MatrixProfile")
    expect_length(result_query, 12)
    expect_true(attr(result_query, "join"))
    expect_equal(result_query$ez, 0)
    expect_equal(round(mean(result_query$mp), 4), 2.8817)
    expect_equal(round(sd(result_query$mp), 4), 0.7918)
    expect_equal(round(mean(result_query$pi), 4), 314.6084)
    expect_equal(round(sd(result_query$pi), 4), 154.1199)
  })

  # multiple windows
  #
  result_multiple <- analyze(ts, windows = c(10:70))

  test_that("Result Multiple", {
    expect_named(result_multiple)
    expect_type(result_multiple, "list")
    expect_s3_class(result_multiple, "PMP")
    expect_length(result_multiple, 9)
    expect_false(attr(result_multiple, "join"))
    expect_equal(result_multiple$ez, 0.5)
    expect_equal(result_multiple$upper_window, 20)
    expect_equal(round(mean(unlist(result_multiple$pmp)), 4), 1.6873)
    expect_equal(round(sd(unlist(result_multiple$pmp)), 3), 0.717)
    expect_equal(round(mean(unlist(result_multiple$pmpi)), 4), 266.6235)
    expect_equal(round(sd(unlist(result_multiple$pmpi)), 4), 151.8338)
  })
}
