context("Testing Annotation functions")
library(tsmp)

data <- as.matrix(mp_test_data$train$data[1:1000])
window <- 50
stop_loc <- 150
prof_size <- nrow(data) - window + 1
comp <- hard <- motion <- stopw <- zero <- NULL

test_that("Error", {
  # window too big
  expect_error(av_complexity(data, nrow(data)), regexp = "too short relative")
  expect_error(av_hardlimit_artifact(data, nrow(data)), regexp = "too short relative")
  expect_error(av_motion_artifact(data, nrow(data)), regexp = "too short relative")
  expect_error(av_stop_word(data, nrow(data), stop_loc), regexp = "too short relative")
  expect_error(av_zerocrossing(data, nrow(data)), regexp = "too short relative")
  # window too small
  expect_error(av_complexity(data, 1), regexp = "be at least 4")
  expect_error(av_hardlimit_artifact(data, 1), regexp = "be at least 4")
  expect_error(av_motion_artifact(data, 1), regexp = "be at least 4")
  expect_error(av_stop_word(data, 1, stop_loc), regexp = "be at least 4")
  expect_error(av_zerocrossing(data, 1), regexp = "be at least 4")
})

test_that("Silent", {
  expect_silent(comp <<- av_complexity(data, window))
  expect_silent(hard <<- av_hardlimit_artifact(data, window))
  expect_silent(motion <<- av_motion_artifact(data, window))
  expect_silent(stopw <<- av_stop_word(data, window, stop_loc))
  expect_silent(zero <<- av_zerocrossing(data, window))
})

test_that("Result dim", {
  expect_equal(dim(comp), c(prof_size, 1))
  expect_equal(dim(hard), c(prof_size, 1))
  expect_equal(dim(motion), c(prof_size, 1))
  expect_equal(dim(stopw), c(prof_size, 1))
  expect_equal(dim(zero), c(prof_size, 1))
})
test_that("Result values", {
  expect_equal(round(sum(comp) / sd(comp), 2), 1689.92)
  expect_equal(round(sum(hard) / sd(hard), 2), 3568.52)
  expect_equal(round(sum(motion) / sd(motion), 1), 1015.7)
  expect_equal(round(sum(stopw) / sd(stopw), 2), 1336.86)
  expect_equal(round(sum(zero) / sd(zero), 2), 666.75)
})
