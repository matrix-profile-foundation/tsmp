if (!testthat:::on_cran()) {
  skip_if_not_installed("RJSONIO")
  library(tsmp)

  result <- compute(mp_toy_data$data[, 1], 80)
  write_json(result, file = file.path(tempdir(), "output.json"))
  result2 <- read_json(file.path(tempdir(), "output.json"))
  unlink(file.path(tempdir(), "output.json"))

  test_that(
    "Reserializing MP",
    expect_true(all.equal(result, result2))
  )

  result <- compute(mp_toy_data$data[, 1])

  write_json(result, file = file.path(tempdir(), "output.json"))
  result2 <- read_json(file.path(tempdir(), "output.json"))
  unlink(file.path(tempdir(), "output.json"))

  test_that(
    "Reserializing PMP",
    expect_true(all.equal(result, result2))
  )
}

test_that("JSON preserves profile vectors and provenance", {
  skip_if_not_installed("RJSONIO")
  file <- tempfile(fileext = ".json")
  on.exit(unlink(file))
  data <- mp_toy_data$data[1:150, 1]

  for (query in list(NULL, mp_toy_data$data[1:100, 2])) {
    result <- compute(data, 30, query = query)
    tsmp::write_json(result, file)
    restored <- tsmp::read_json(file)
    expect_equal(restored, result)
    expect_null(dim(restored$mp))
    expect_null(dim(restored$pi))
    expect_equal(attr(restored, "origin"), attr(result, "origin"))
    expect_identical(attr(restored, "partial"), attr(result, "partial"))
  }
})

test_that("JSON files without provenance remain readable", {
  skip_if_not_installed("RJSONIO")
  file <- tempfile(fileext = ".json")
  on.exit(unlink(file))
  result <- compute(mp_toy_data$data[1:150, 1], 30)
  attr(result, "origin") <- NULL
  attr(result, "partial") <- NULL
  tsmp::write_json(result, file)
  expect_equal(tsmp::read_json(file), result)
})

test_that("write_json rejects unsupported objects", {
  expect_error(tsmp::write_json("text", tempfile()), "MatrixProfile or PMP")
})

test_that("failed JSON writes restore the digits option", {
  skip_if_not_installed("RJSONIO")
  result <- compute(mp_toy_data$data[1:150, 1], 30)
  digits <- getOption("digits")
  expect_error(suppressWarnings(tsmp::write_json(result, file.path(tempfile(), "missing.json"))))
  expect_identical(getOption("digits"), digits)
})
