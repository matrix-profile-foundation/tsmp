context("Testing Annotation functions")
library(tsmp)

data <- as.matrix(test_data$train$data[1:1000])
window <- 50
stop.loc <- 150
prof.size <- nrow(data) - window + 1
comp <- hard <- motion <- stopw <- zero <- NULL

test_that("Error", {
  # window too big
  expect_error(av.complexity(data, nrow(data)), regexp = "too short relative")
  expect_error(av.hardlimit.artifact(data, nrow(data)), regexp = "too short relative")
  expect_error(av.motion.artifact(data, nrow(data)), regexp = "too short relative")
  expect_error(av.stop.word(data, nrow(data), stop.loc), regexp = "too short relative")
  expect_error(av.zerocrossing(data, nrow(data)), regexp = "too short relative")
  # window too small
  expect_error(av.complexity(data, 1), regexp = "be at least 4")
  expect_error(av.hardlimit.artifact(data, 1), regexp = "be at least 4")
  expect_error(av.motion.artifact(data, 1), regexp = "be at least 4")
  expect_error(av.stop.word(data, 1, stop.loc), regexp = "be at least 4")
  expect_error(av.zerocrossing(data, 1), regexp = "be at least 4")
})

test_that("Silent", {
  expect_silent(comp <<- av.complexity(data, window))
  expect_silent(hard <<- av.hardlimit.artifact(data, window))
  expect_silent(motion <<- av.motion.artifact(data, window))
  expect_silent(stopw <<- av.stop.word(data, window, stop.loc))
  expect_silent(zero <<- av.zerocrossing(data, window))
})

test_that("Result dim", {
  expect_equal(dim(comp), c(prof.size, 1))
  expect_equal(dim(hard), c(prof.size, 1))
  expect_equal(dim(motion), c(prof.size, 1))
  expect_equal(dim(stopw), c(prof.size, 1))
  expect_equal(dim(zero), c(prof.size, 1))
})
test_that("Result values", {
  expect_equal(round(sum(comp) / sd(comp), 2), 1689.92)
  expect_equal(round(sum(hard) / sd(hard), 2), 3568.52)
  expect_equal(round(sum(motion) / sd(motion), 1), 1015.7)
  expect_equal(round(sum(stopw) / sd(stopw), 2), 1336.86)
  expect_equal(round(sum(zero) / sd(zero), 2), 666.75)
})
