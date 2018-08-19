context("Testing FLUSS functions")
library(tsmp)

data <- fluss_data$tilt.abp$data[1:1000]
w <- 10
truth <- 400
nseg <- 3
mp <- mstomp(data, w, verbose = 0)
cac <- fluss.cac(mp$pi, w)
segments <- fluss.extract(cac, nseg, w)
score <- fluss.score(truth, segments, length(data))
res <- fluss(t(data), w, nseg, gtruth = truth, verbose = 0)
res.nt <- fluss(data, w, nseg, verbose = 0)

test_that("Errors", {
  # big window size
  expect_error(fluss(table(data), w, nseg, gtruth = truth, verbose = 0), regexp = "Unknown type of data")
})

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

test_that("Full fluss", {
  expect_equal(res$score, score)
  expect_equal(res$segments, segments)
  expect_equal(res$segments, res.nt$segments)
  expect_equal(res$mp, res.nt$mp)
  expect_equal(res$pi, res.nt$pi)
  expect_equal(res$cac, res.nt$cac)
  expect_equal(round(mean(res$cac), 4), round(mean(cac), 4))
  expect_equal(round(mean(res$cac), 4), round(mean(res.nt$cac), 4))
  expect_equal(round(sd(res$cac), 4), round(sd(cac), 4))
  expect_equal(round(sd(res$cac), 4), round(sd(res.nt$cac), 4))
  expect_equal(round(min(res$cac), 4), round(min(cac), 4))
  expect_equal(round(min(res$cac), 4), round(min(res.nt$cac), 4))
  expect_equal(max(res$cac), max(cac))
  expect_equal(max(res$cac), max(res.nt$cac))
})
