if (!testthat:::on_cran()) {
  context("Testing Serialization")
  library(tsmp)

  result <- compute(mp_toy_data$data[, 1], 80)
  write(result, file = file.path(tempdir(), "output.json"))
  result2 <- read(file.path(tempdir(), "output.json"))
  unlink(file.path(tempdir(), "output.json"))

  test_that(
    "Reserializing MP",
    expect_true(all.equal(result, result2))
  )

  result <- compute(mp_toy_data$data[, 1])

  write(result, file = file.path(tempdir(), "output.json"))
  result2 <- read(file.path(tempdir(), "output.json"))
  unlink(file.path(tempdir(), "output.json"))

  test_that(
    "Reserializing PMP",
    expect_true(all.equal(result, result2))
  )
}
