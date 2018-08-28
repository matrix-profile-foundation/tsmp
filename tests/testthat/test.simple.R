context("Testing SiMPle Fast")
library(tsmp)

w <- 30
data <- toy_data$data[1:250, ] # 3 dimensions matrix
query <- toy_data$data[251:500, ] # 3 dimensions matrix

test_that("Errors", {
  # big window size
  expect_error(simple.fast(data, window.size = 500), regexp = "too short relative")

  # short window size
  expect_error(simple.fast(data, window.size = 2), regexp = "must be at least 4")

  # invalid window
  expect_error(simple.fast(data, window.size = data), regexp = "type of window.size")

  # data and query dim must be the same
  expect_error(simple.fast(data, data[, 1], window.size = w), regexp = "Data and query dimensions")
  expect_error(simple.fast(data[, 1], data, window.size = w), regexp = "Data and query dimensions")

  # Unknown type
  expect_error(simple.fast(table(data), window.size = w), regexp = "Unknown type of data")
  expect_error(simple.fast(data, table(data), window.size = w), regexp = "Unknown type of query")
})

test_that("Messages", {
  expect_message(simple.fast(data, window.size = w, verbose = 1), regexp = "Finished")
})

if (skip_on_cran()) {
  result.self <- simple.fast(list(data[, 1], data[, 2], data[, 3]), window.size = w, verbose = 2)
  result.join <- simple.fast(as.data.frame(t(data)), as.data.frame(t(query)), window.size = w, verbose = 2)
} else {
  result.self <- simple.fast(data, window.size = w, verbose = 0)
  result.join <- simple.fast(data, query, window.size = w, verbose = 0)
}

test_that("SiMPle Results", {
  expect_equal(round(sum(result.self$mp), 3), 419.509)
  expect_equal(round(sd(result.self$mp), 3), 0.841)
  expect_equal(sum(result.self$pi), 23878)
  expect_equal(round(sd(result.self$pi), 3), 64.977)
  expect_equal(round(sum(result.join$mp), 3), 908.248)
  expect_equal(round(sd(result.join$mp), 3), 2.263)
  expect_equal(sum(result.join$pi), 24981)
  expect_equal(round(sd(result.join$pi), 3), 61.021)
})
