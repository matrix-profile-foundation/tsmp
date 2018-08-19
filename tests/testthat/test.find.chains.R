context("Testing Time Series Chains")
library(tsmp)

w <- 50
data <- gait_data
mp <- stamp(data, window.size = w, exclusion.zone = 1 / 4, verbose = 0)
res <- find.chains(mp)

test_that("Find Chains", {
  expect_equal(length(res), 2)
  expect_equal(length(res$chains), 58)
  expect_equal(length(res$best), 6)
  expect_known_hash(res$chains, "d7c3f43152")
})
