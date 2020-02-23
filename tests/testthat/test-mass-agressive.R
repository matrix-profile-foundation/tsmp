if (!testthat:::on_cran()) {
  context("Testing if mass v2 and mass v3 agrees")
  library(tsmp)

  ref_data <- mp_fluss_data$tilt_abp$data[1:1000]
  # minimum example, data and query

  windows <- sample((length(ref_data) %/% 3), 100)
  windows <- windows[windows > 4]

  for (w in windows) {
    for (i in sample(length(ref_data), 10)) {
      nn2 <- dist_profile(ref_data, ref_data, window_size = w, index = i, method = "v2")
      nn3 <- dist_profile(ref_data, ref_data, window_size = w, index = i, method = "v3")
      alleq <- all(
        all.equal(nn2$distance_profile, nn3$distance_profile),
        all.equal(nn2$last_product, nn3$last_product),
        all.equal(nn2$par$data_fft, nn3$par$data_fft)
      )

      test_that("Final result", {
        expect_true(alleq)
      })
    }
  }
}
