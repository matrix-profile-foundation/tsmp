if (!testthat:::on_cran()) {
  context("Testing Serialization")
  library(tsmp)

  result <- compute(mp_toy_data$data[, 1], 80)
  write(result, file = "output.json")
  result2 <- result <- read("output.json")
  unlink("output.json")

  test_that(
    "reserializing",
    expect_true(all.equal(result, result2))
  )
}
