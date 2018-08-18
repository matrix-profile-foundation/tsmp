context("Testing FLUSS functions")
library(tsmp)

data <- fluss_data$tilt.abp$data[1:1000]
w <- 10
truth <- 400
mp <- mstomp(data, w, verbose = 0)
cac <- fluss.cac(mp$pi, w)
segments <- fluss.extract(cac, 3, w)
score <- fluss.score(truth, segments, length(data))

test_that("Corrected Arc Count", {
  expect_equal(round(mean(cac), 4), 0.9941)
  expect_equal(round(sd(cac), 4), 0.0187)
  expect_equal(round(min(cac), 4), 0.8838)
  expect_equal(max(cac), 1)
})

test_that("Segments found", {
  expect_equal(segments, c(941, 875, 141))
})

test_that("Score", {
  expect_equal(score, 0.259)
})
