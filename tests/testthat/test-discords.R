if (!testthat:::on_cran()) {
  context("Testing API discords")
  library(tsmp)

  ts <- mp_toy_data$data[, 1]

  result <- compute(ts, windows = 30) %>% discords()

  test_that("Discords", {
    expect_named(result)
    expect_type(result, "list")
    expect_s3_class(result, "MatrixProfile")
    expect_length(result, 11)
    expect_false(attr(result, "join"))
    expect_equal(round(mean(unlist(result$discord)), 4), 298.8182)
    expect_equal(round(sd(unlist(result$discord)), 4), 148.6915)
  })
}
