context("Testing Annotation functions")
library(tsmp)

data <- finger3500
window <- 3500
prof.size <- nrow(data) - window + 1
comp <- hard <- motion <- stopw <- zero <- NULL

test_that("Error", {
  # window too big
  expect_error(av.complexity(data, nrow(data)), regexp = "too short relative")
  expect_error(av.hardlimit.artifact(data, nrow(data)), regexp = "too short relative")
  expect_error(av.motion.artifact(data, nrow(data)), regexp = "too short relative")
  expect_error(av.stop.word(data, nrow(data)), regexp = "too short relative")
  expect_error(av.zerocrossing(data, nrow(data)), regexp = "too short relative")
  # window too small
  expect_error(av.complexity(data, 1), regexp = "be at least 4")
  expect_error(av.hardlimit.artifact(data, 1), regexp = "be at least 4")
  expect_error(av.motion.artifact(data, 1), regexp = "be at least 4")
  expect_error(av.stop.word(data, 1), regexp = "be at least 4")
  expect_error(av.zerocrossing(data, 1), regexp = "be at least 4")
})

test_that("Silent", {
  expect_silent(comp <<- av.complexity(data, window))
  expect_silent(hard <<- av.hardlimit.artifact(data, window))
  expect_silent(motion <<- av.motion.artifact(data, window))
  expect_silent(stopw <<- av.stop.word(data, window))
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
  expect_equal(round(sum(comp) / sd(comp), 2), 40637.64)
  expect_equal(round(sum(hard) / sd(hard), 1), 150580.2)
  expect_equal(round(sum(motion) / sd(motion), 2), 80083.71)
  expect_equal(round(sum(stopw) / sd(stopw), 2), 41982.67)
  expect_equal(round(sum(zero) / sd(zero), 2), 42503.32)
})
