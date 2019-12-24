if (!testthat:::on_cran()) {

  # Stamps and Stomps agree ----

  context("Testing if Stamps and Stomps algorithms agree")
  library(tsmp)

  ## Test errors ----

  test_that("Errors", {
    # big window size
    expect_error(mstomp(mp_toy_data$data[1:400, ], window_size = 500), "too short relative")
    expect_error(mstomp_par(mp_toy_data$data[1:400, ], window_size = 500), "too short relative")
    expect_error(stomp(mp_toy_data$data[1:400, ], window_size = 500), "too short relative")
    expect_error(stomp_par(mp_toy_data$data[1:400, ], window_size = 500), "too short relative")
    expect_error(stamp(mp_toy_data$data[1:400, ], window_size = 500), "too short relative")
    expect_error(stamp_par(mp_toy_data$data[1:400, ], window_size = 500), "too short relative")
    expect_error(scrimp(mp_toy_data$data[1:400, ], window_size = 500), "too short relative")

    # intersect
    expect_error(mstomp(mp_toy_data$data[1:400, ], window_size = 30, must_dim = c(1, 2), exc_dim = c(2, 3)), "presented in both")
    expect_error(mstomp_par(mp_toy_data$data[1:400, ], window_size = 30, must_dim = c(1, 2), exc_dim = c(2, 3)), "presented in both")
    # too many must_dim
    expect_error(mstomp(mp_toy_data$data[1:400, ], window_size = 30, must_dim = c(1, 2, 3, 4)), "must_dim")
    expect_error(mstomp_par(mp_toy_data$data[1:400, ], window_size = 30, must_dim = c(1, 2, 3, 4)), "must_dim")
    # too many exc_dim
    expect_error(mstomp(mp_toy_data$data[1:400, ], window_size = 30, exc_dim = c(1, 2, 3, 4)), "exc_dim")
    expect_error(mstomp_par(mp_toy_data$data[1:400, ], window_size = 30, exc_dim = c(1, 2, 3, 4)), "exc_dim")

    # small window size
    expect_error(stamp(mp_toy_data$data[1:400, ], window_size = 2), "window_size")
    expect_error(stamp_par(mp_toy_data$data[1:400, ], window_size = 2), "window_size")
    expect_error(mstomp(mp_toy_data$data[1:400, ], window_size = 2), "window_size")
    expect_error(mstomp_par(mp_toy_data$data[1:400, ], window_size = 2), "window_size")
    expect_error(stomp(mp_toy_data$data[1:400, ], window_size = 2), "window_size")
    expect_error(stomp_par(mp_toy_data$data[1:400, ], window_size = 2), "window_size")
    expect_error(scrimp(mp_toy_data$data[1:400, ], window_size = 2), "window_size")

    # unknown data type
    expect_error(stamp(table(rpois(100, 5)), window_size = 30), "Unknown type")
    expect_error(stamp_par(table(rpois(100, 5)), window_size = 30), "Unknown type")
    expect_error(mstomp(table(rpois(100, 5)), window_size = 30), "Unknown type")
    expect_error(mstomp_par(table(rpois(100, 5)), window_size = 30), "Unknown type")
    expect_error(stomp(table(rpois(100, 5)), window_size = 30), "Unknown type")
    expect_error(stomp_par(table(rpois(100, 5)), window_size = 30), "Unknown type")
    expect_error(scrimp(table(rpois(100, 5)), window_size = 30), "Unknown type")
  })

  ## Test finish ----

  test_that("Finish", {
    expect_message(tsmp(mp_toy_data$data[1:400, 1], mode = "stamp", window_size = 30), "Finished")
    expect_message(tsmp(mp_toy_data$data[1:400, 1], mode = "stamp", window_size = 30, n_workers = 2), "Finished")
    expect_message(tsmp(mp_toy_data$data[1:400, 1], mode = "stomp", window_size = 30), "Finished")
    expect_message(tsmp(mp_toy_data$data[1:400, 1], mode = "stomp", window_size = 30, n_workers = 2), "Finished")
    expect_message(tsmp(mp_toy_data$data[1:400, 1], mode = "mstomp", window_size = 30), "Finished")
    expect_message(tsmp(mp_toy_data$data[1:400, 1], mode = "mstomp", window_size = 30, n_workers = 2), "Finished")
    expect_message(tsmp(mp_toy_data$data[1:400, 1], mode = "scrimp", window_size = 30), "Finished")
    expect_message(tsmp(mp_toy_data$data[1:400, 1], mp_toy_data$data[1:400, 1], mode = "scrimp", window_size = 30), "not implemented")
  })

  ## Create MP's ----

  # STAMP
  stamp_test <- stamp(mp_toy_data$data[1:400, 1], window_size = 30, verbose = 0)
  stamp_join_test <- stamp(mp_toy_data$data[1:400, 1], mp_toy_data$data[1:100, 2], window_size = 30, verbose = 0)
  stamp_par_test <- stamp_par(mp_toy_data$data[1:400, 1], window_size = 30, verbose = 0)
  stamp_par_join_test <- stamp_par(mp_toy_data$data[1:400, 1], mp_toy_data$data[1:100, 2], window_size = 30, verbose = 0)

  # STOMP
  stomp_test <- stomp(mp_toy_data$data[1:400, 1], window_size = 30, verbose = 0)
  stompi_test <- tsmp(mp_toy_data$data[1:200, 1], window_size = 30, verbose = 0)
  stompi_test <- stompi_update(stompi_test, mp_toy_data$data[201:400, 1])
  stomp_join_test <- stomp(mp_toy_data$data[, 1], mp_toy_data$data[1:400, 2], window_size = 30, verbose = 0)
  stomp_par_test <- stomp_par(mp_toy_data$data[1:400, 1], window_size = 30, verbose = 0)
  stomp_par_join_test <- stomp_par(mp_toy_data$data[, 1], mp_toy_data$data[1:400, 2], window_size = 30, verbose = 0)

  # MSTOMP Uni
  mstomp_test1 <- mstomp(mp_toy_data$data[1:400, 1], window_size = 30, verbose = 0)
  mstomp_par_test1 <- mstomp_par(mp_toy_data$data[1:400, 1], window_size = 30, verbose = 0)
  # MSTOMP Multi
  mstomp_test <- mstomp(mp_toy_data$data[1:400, ], window_size = 30, verbose = 0)
  mstomp_test_must <- mstomp(mp_toy_data$data[1:400, ], window_size = 30, must_dim = c(1, 2), verbose = 0)
  mstomp_test_exc <- mstomp(mp_toy_data$data[1:400, ], window_size = 30, exc_dim = c(1, 2), verbose = 0)
  mstomp_par_test <- mstomp_par(mp_toy_data$data[1:400, ], window_size = 30, verbose = 0)
  mstomp_par_test_must <- mstomp_par(mp_toy_data$data[1:400, ], window_size = 30, must_dim = c(1, 2), verbose = 0)
  mstomp_par_test_exc <- mstomp_par(mp_toy_data$data[1:400, ], window_size = 30, exc_dim = c(1, 2), verbose = 0)

  scrimp_test <- scrimp(mp_toy_data$data[1:400, 1], window_size = 30, verbose = 0)

  ## Check consistency ----

  test_that("Basic Results", {
    expect_equal(round(sum(stamp_test$mp) / sd(stamp_test$mp), 3), 1117.483)
    expect_equal(sum(which(is.infinite(stamp_test$rmp))), 5816)
    expect_equal(round(sum(stamp_test$rmp[1:155]) / sd(stamp_test$rmp[1:155]), 3), 451.847)
    expect_equal(sum(which(is.infinite(stamp_test$lmp))), 136)
    expect_equal(round(sum(stamp_test$lmp[17:171]) / sd(stamp_test$lmp[17:171]), 3), 345.977)
    expect_equal(round(sum(stamp_test$pi) / sd(stamp_test$pi), 3), 567.591)
    expect_equal(round(sum(stamp_test$rpi[1:355]) / sd(stamp_test$rpi[1:355]), 3), 1354.43)
    expect_equal(round(sum(stamp_test$lpi[17:371]) / sd(stamp_test$lpi[17:371]), 3), 340.923)
    expect_equal(stamp_test$w, 30)
    expect_equal(stamp_test$ez, 0.5)
    expect_equal(class(stamp_test), "MatrixProfile")
    expect_equal(class(stomp_test), "MatrixProfile")
    expect_equal(class(mstomp_test), "MultiMatrixProfile")
  })

  ## Check STOMP and STOMPI ----

  test_that("Stompi Results", {
    expect_equal(stomp_test$mp, stompi_test$mp)
    expect_equal(stomp_test$pi, stompi_test$pi)
    expect_equal(stomp_test$rmp, stompi_test$rmp)
    expect_equal(stomp_test$rpi, stompi_test$rpi)
    expect_equal(stomp_test$lmp, stompi_test$lmp)
    expect_equal(stomp_test$lpi, stompi_test$lpi)
  })

  ## Check consistency for SCRIMP ----

  test_that("Scrimp Results", {
    expect_equal(class(scrimp_test), "MatrixProfile")
    expect_equal(round(sum(scrimp_test$mp) / sd(scrimp_test$mp), 2), 1117.48)
    expect_equal(round(sum(scrimp_test$pi) / sd(scrimp_test$pi), 3), 567.591)
    expect_equal(scrimp_test$w, 30)
    expect_equal(scrimp_test$ez, 0.5)
  })

  ## Check inter-consistency ----

  # stamp_test and stamp_par_test
  test_that("Stamp equals to Stamp_par", {
    expect_equal(stamp_test, stamp_par_test)
  })

  # stamp_join_test and stamp_par_join_test
  test_that("Stamp Join equals to Stamp_par Join", {
    expect_equal(stamp_join_test, stamp_par_join_test)
  })

  # stamp_test and stomp_test
  test_that("Stamp equals to Stomp", {
    expect_equal(stamp_test, stomp_test)
  })

  # scrimp_test and stomp_test
  test_that("Scrimp equals to Stomp", {
    expect_equal(scrimp_test, stomp_test)
  })

  test_that("Stomp Join equals to Stomp_par Join", {
    expect_equal(stomp_join_test, stomp_par_join_test)
  })

  # stamp_test and stomp_par_test
  test_that("Stamp equals to Stomp_par", {
    expect_equal(stamp_test, stomp_par_test)
  })

  # stamp_test and mstomp_test1
  test_that("Stamp equals to mStomp", {
    expect_equal(stamp_test, mstomp_test1)
  })

  # stamp_test and mstomp_par_test1
  test_that("Stamp equals to Stomp", {
    expect_equal(stamp_test, mstomp_par_test1)
  })

  # mstomp_test and mstomp_par_test
  test_that("mStomp equals to mStomp_par", {
    expect_equal(mstomp_test, mstomp_par_test)
  })

  # mstomp_test_must and mstomp_par_test_must
  test_that("mStomp must equals to mStomp_par must", {
    expect_equal(mstomp_test_must, mstomp_par_test_must)
  })

  # mstomp_test_exc and mstomp_par_test_exc
  test_that("mStomp exc equals to mStomp_par exc", {
    expect_equal(mstomp_test_exc, mstomp_par_test_exc)
  })
}
