if (!testthat:::on_cran()) {
  context("Testing API motifs")
  library(tsmp)

  ts <- mp_toy_data$data[, 1]

  result <- compute(ts, windows = 30) %>% motifs()

  test_that("Motifs", {
    expect_named(result)
    expect_type(result, "list")
    expect_s3_class(result, "MatrixProfile")
    expect_length(result, 11)
    expect_false(attr(result, "join"))
    expect_equal(round(mean(unlist(result$motif)), 4), 256.3182)
    expect_equal(round(sd(unlist(result$motif)), 4), 167.7008)
  })
}
