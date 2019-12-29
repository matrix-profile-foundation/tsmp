if (!testthat:::on_cran()) {
  context("Testing MPdist")
  library(tsmp)

  ref_data <- mp_toy_data$data[, 1]
  qe_data <- mp_toy_data$data[, 2]
  qd_data <- mp_toy_data$data[150:200, 1]
  w <- mp_toy_data$sub_len

  # distance between data of same size
  deq <- mpdist(ref_data, qe_data, w)

  # distance between data of different sizes
  ddiff <- mpdist(qe_data, qd_data, w)

  # distance vector between data of different sizes
  ddvect <- mpdist(ref_data, qd_data, w, type = "vector")

  test_that("MPdist", {
    expect_equal(round(deq, 5), 2.02497)
    expect_equal(round(ddiff, 5), 5.69151)
    expect_equal(round(mean(ddvect$mpdist), 5), 4.63318)
    expect_equal(round(sd(ddvect$mpdist), 5), 1.55205)
  })
}
