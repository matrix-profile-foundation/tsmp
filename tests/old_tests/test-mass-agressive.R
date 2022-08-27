


test_that("MASS stress test", {
  ref_data <- as.numeric(head(mp_fluss_data$tilt_abp$data, 1000))
  # minimum example, data and query

  set.seed(2022)
  windows <- sample((length(ref_data) %/% 3), 100)
  windows <- windows[windows > 4]

  set.seed(2200)
  for (w in windows) {
    for (i in sample(length(ref_data), 10)) {
      nn2 <- dist_profile(ref_data, ref_data, window_size = w, index = i)
      nn3 <- dist_profile(ref_data, ref_data, window_size = w, index = i)
      alleq <- all(
        all.equal(nn2$distance_profile, nn3$distance_profile),
        all.equal(nn2$last_product, nn3$last_product),
        all.equal(nn2$par$data_fft, nn3$par$data_fft)
      )

      expect_true(alleq)
    }
  }
})
