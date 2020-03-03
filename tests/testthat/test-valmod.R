if (!testthat:::on_cran()) {
  context("Testing Valmod function")

  ref_data <- tsmp::mp_toy_data$data[, 1L]
  query_data <- tsmp::mp_toy_data$data[, 2L]
  # self similarity
  mp <- valmod(ref_data, window_min = 30L, window_max = 40L)
  # join similarity
  mp_join <- valmod(ref_data, query_data, window_min = 30L, window_max = 40L, verbose = 1L)

  test_that("Success", {
    expect_silent(valmod(ref_data, window_min = 30L, window_max = 40L, verbose = 0L))
    expect_message(valmod(ref_data, window_min = 30L, window_max = 40L, verbose = 1L))
  })

  test_that("Self similarity", {
    expect_equal(round(sum(mp$mp) / sd(mp$mp), 3L), 1584.183)
    expect_equal(round(sum(mp$pi) / sd(mp$pi), 3L), 741.849)
    expect_equal(round(sum(mp$w) / sd(mp$w), 2L), 11406.86)
  })

  test_that("Join similarity", {
    expect_equal(round(sum(mp_join$mp) / sd(mp_join$mp), 3L), 2190.417)
    expect_equal(round(sum(mp_join$pi) / sd(mp_join$pi), 2L), 922.14)
    expect_equal(round(sum(mp_join$w) / sd(mp_join$w), 1L), 10532.1)
  })
}
