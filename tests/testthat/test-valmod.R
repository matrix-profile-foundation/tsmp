if (skip_on_cran()) {
  context("Testing Valmod function")

  ref_data <- tsmp::mp_toy_data$data[, 1]
  query_data <- tsmp::mp_toy_data$data[, 2]
  # self similarity
  mp <- valmod(ref_data, window_min = 30, window_max = 40)
  # join similarity
  mp_join <- valmod(ref_data, query_data, window_min = 30, window_max = 40, verbose = 1)

  test_that("Success", {
    expect_silent(valmod(ref_data, window_min = 30, window_max = 40, verbose = 0))
    expect_message(valmod(ref_data, window_min = 30, window_max = 40, verbose = 1))
  })

  test_that("Self similarity", {
    expect_equal(round(sum(mp$mp) / sd(mp$mp), 3), 1584.183)
    expect_equal(round(sum(mp$pi) / sd(mp$pi), 3), 741.849)
    expect_equal(round(sum(mp$w) / sd(mp$w), 2), 11406.86)
  })

  test_that("Join similarity", {
    expect_equal(round(sum(mp_join$mp) / sd(mp_join$mp), 3), 2190.417)
    expect_equal(round(sum(mp_join$pi) / sd(mp_join$pi), 2), 922.14)
    expect_equal(round(sum(mp_join$w) / sd(mp_join$w), 1), 10532.1)
  })
}
