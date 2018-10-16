if (identical(Sys.getenv("STRESS_TEST"), "true") &&
  identical(Sys.getenv("NOT_CRAN"), "true") &&
  !identical(Sys.getenv("TRAVIS"), "true") &&
  !identical(Sys.getenv("APPVEYOR"), "True")) {
  context("Stress test for STAMP, STOMP and SCRIMP")
  library(tsmp)

  ssize <- 10

  set.seed(2018)
  windows <- c(4, sample(50:4000, (ssize - 1)))
  message("Windows: ", windows)

  for (w in windows) {
    min_data_size <- 2 * w + 1
    max_data_size <- 10000

    data_sizes <- sample(seq(min_data_size, max_data_size), ssize)
    message("Window: ", w)
    message("Sizes: ", data_sizes)

    for (ds in data_sizes) {
      message("Window: ", w)
      message("Size: ", ds)

      data1 <- mp_fluss_data$tilt_abp$data[1:ds] # 40000
      data2 <- mp_fluss_data$walkjogrun$data[1:ds] # 10000
      data3 <- mp_meat_data$sub$data[1:ds] # 107520 with RW
      data4 <- mp_test_data$train$data[1:ds] # 215010

      stamp_test_1 <- stamp_par(data1, window_size = w, n_workers = 4)
      stamp_test_2 <- stamp_par(data2, window_size = w, n_workers = 4)
      stamp_test_3 <- stamp_par(data3, window_size = w, n_workers = 4)
      stamp_test_4 <- stamp_par(data4, window_size = w, n_workers = 4)

      stomp_test_1 <- stomp_par(data1, window_size = w, n_workers = 4)
      stomp_test_2 <- stomp_par(data2, window_size = w, n_workers = 4)
      stomp_test_3 <- stomp_par(data3, window_size = w, n_workers = 4)
      stomp_test_4 <- stomp_par(data4, window_size = w, n_workers = 4)

      scrimp_test_1 <- scrimp(data1, window_size = w)
      scrimp_test_2 <- scrimp(data2, window_size = w)
      scrimp_test_3 <- scrimp(data3, window_size = w)
      scrimp_test_4 <- scrimp(data4, window_size = w)

      label <- paste("Window", w, "data size", ds)
      test_that(label, {
        expect_true(all.equal(stamp_test_1$mp, stomp_test_1$mp, tolerance = 0.01), info = paste("stamp, stomp, 1", label))
        expect_true(all.equal(stamp_test_2$mp, stomp_test_2$mp, tolerance = 0.01), info = paste("stamp, stomp, 2", label))
        expect_true(all.equal(stamp_test_3$mp, stomp_test_3$mp, tolerance = 0.01), info = paste("stamp, stomp, 3", label))
        expect_true(all.equal(stamp_test_4$mp, stomp_test_4$mp, tolerance = 0.01), info = paste("stamp, stomp, 4", label))

        expect_true(all.equal(stamp_test_1$mp, scrimp_test_1$mp, tolerance = 0.01), info = paste("stamp, scrimp, 1", label))
        expect_true(all.equal(stamp_test_2$mp, scrimp_test_2$mp, tolerance = 0.01), info = paste("stamp, scrimp, 2", label))
        expect_true(all.equal(stamp_test_3$mp, scrimp_test_3$mp, tolerance = 0.01), info = paste("stamp, scrimp, 3", label))
        expect_true(all.equal(stamp_test_4$mp, scrimp_test_4$mp, tolerance = 0.01), info = paste("stamp, scrimp, 4", label))

        expect_true(all.equal(stomp_test_1$mp, scrimp_test_1$mp, tolerance = 0.01), info = paste("stomp, scrimp, 1", label))
        expect_true(all.equal(stomp_test_2$mp, scrimp_test_2$mp, tolerance = 0.01), info = paste("stomp, scrimp, 2", label))
        expect_true(all.equal(stomp_test_3$mp, scrimp_test_3$mp, tolerance = 0.01), info = paste("stomp, scrimp, 3", label))
        expect_true(all.equal(stomp_test_4$mp, scrimp_test_4$mp, tolerance = 0.01), info = paste("stomp, scrimp, 4", label))
      })
    }
  }
}
