context("Testing SiMPle Fast")
library(tsmp)

w <- 30
data <- toy_data$data # 3 dimensions matrix
if (skip_on_cran()) {
  result <- simple.fast(data, w, verbose = 2)
} else {
  result <- simple.fast(data, w, verbose = 0)
}

test_that("SiMPle Results", {
  expect_equal(round(sum(result$mp), 3), 806.132)
  expect_equal(round(sd(result$mp), 3), 0.575)
  expect_equal(sum(result$pi), 135450)
  expect_equal(round(sd(result$pi), 3), 151.06)
})
