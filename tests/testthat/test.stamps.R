context("Testing if Stamps and Stomps algorithms agree")
library(tsmp)

if (skip_on_cran()) {
  test_that("Errors", {
    # big window size
    expect_error(mstomp(mp_toy_data$data[1:200, ], window_size = 500), regexp = "too short relative")
    expect_error(mstomp_par(mp_toy_data$data[1:200, ], window_size = 500), regexp = "too short relative")
    expect_error(stomp(mp_toy_data$data[1:200, ], window_size = 500), regexp = "too short relative")
    expect_error(stomp_par(mp_toy_data$data[1:200, ], window_size = 500), regexp = "too short relative")
    expect_error(stamp(mp_toy_data$data[1:200, ], window_size = 500), regexp = "too short relative")
    expect_error(stamp_par(mp_toy_data$data[1:200, ], window_size = 500), regexp = "too short relative")

    # intersect
    expect_error(mstomp(mp_toy_data$data[1:200, ], window_size = 30, must_dim = c(1, 2), exc_dim = c(2, 3)), regexp = "presented in both")
    expect_error(mstomp_par(mp_toy_data$data[1:200, ], window_size = 30, must_dim = c(1, 2), exc_dim = c(2, 3)), regexp = "presented in both")
    # too many must_dim
    expect_error(mstomp(mp_toy_data$data[1:200, ], window_size = 30, must_dim = c(1, 2, 3, 4)), regexp = "must have dimension must be less")
    expect_error(mstomp_par(mp_toy_data$data[1:200, ], window_size = 30, must_dim = c(1, 2, 3, 4)), regexp = "must have dimension must be less")
    # too many exc_dim
    expect_error(mstomp(mp_toy_data$data[1:200, ], window_size = 30, exc_dim = c(1, 2, 3, 4)), regexp = "exclusion dimension must be less")
    expect_error(mstomp_par(mp_toy_data$data[1:200, ], window_size = 30, exc_dim = c(1, 2, 3, 4)), regexp = "exclusion dimension must be less")

    # small window size
    expect_error(stamp(mp_toy_data$data[1:200, ], window_size = 2), regexp = "Window size")
    expect_error(stamp_par(mp_toy_data$data[1:200, ], window_size = 2), regexp = "Window size")
    expect_error(mstomp(mp_toy_data$data[1:200, ], window_size = 2), regexp = "Window size")
    expect_error(mstomp_par(mp_toy_data$data[1:200, ], window_size = 2), regexp = "Window size")
    expect_error(stomp(mp_toy_data$data[1:200, ], window_size = 2), regexp = "Window size")
    expect_error(stomp_par(mp_toy_data$data[1:200, ], window_size = 2), regexp = "Window size")

    # unknown data type
    expect_error(stamp(table(rpois(100, 5)), window_size = 30), regexp = "Unknown type")
    expect_error(stamp_par(table(rpois(100, 5)), window_size = 30), regexp = "Unknown type")
    expect_error(mstomp(table(rpois(100, 5)), window_size = 30), regexp = "Unknown type")
    expect_error(mstomp_par(table(rpois(100, 5)), window_size = 30), regexp = "Unknown type")
    expect_error(stomp(table(rpois(100, 5)), window_size = 30), regexp = "Unknown type")
    expect_error(stomp_par(table(rpois(100, 5)), window_size = 30), regexp = "Unknown type")
  })

  test_that("Finish", {
    expect_message(stamp(mp_toy_data$data[1:200, 1], window_size = 30), regex = "Finished")
    expect_message(stamp_par(mp_toy_data$data[1:200, 1], window_size = 30), regex = "Finished")
    expect_message(mstomp(mp_toy_data$data[1:200, 1], window_size = 30), regex = "Finished")
    expect_message(mstomp_par(mp_toy_data$data[1:200, 1], window_size = 30), regex = "Finished")
    expect_message(stomp(mp_toy_data$data[1:200, 1], window_size = 30), regex = "Finished")
    expect_message(stomp_par(mp_toy_data$data[1:200, 1], window_size = 30), regex = "Finished")
  })

  # STAMP
  stamp_test <- stamp(mp_toy_data$data[1:200, 1], window_size = 30, verbose = 0)
  stamp_join_test <- stamp(mp_toy_data$data[1:200, 1], mp_toy_data$data[1:100, 2], window_size = 30, verbose = 0)
  stamp_par_test <- stamp_par(mp_toy_data$data[1:200, 1], window_size = 30, verbose = 0)
  stamp_par_join_test <- stamp_par(mp_toy_data$data[1:200, 1], mp_toy_data$data[1:100, 2], window_size = 30, verbose = 0)

  # STOMP
  stomp_test <- stomp(mp_toy_data$data[1:200, 1], window_size = 30, verbose = 0)
  stomp_join_test <- stomp(mp_toy_data$data[1:200, 1], mp_toy_data$data[1:100, 2], window_size = 30, verbose = 0)
  stomp_par_test <- stomp_par(mp_toy_data$data[1:200, 1], window_size = 30, verbose = 0)
  stomp_par_join_test <- stomp_par(mp_toy_data$data[1:200, 1], mp_toy_data$data[1:100, 2], window_size = 30, verbose = 0)

  # MSTOMP Uni
  mstomp_test1 <- mstomp(mp_toy_data$data[1:200, 1], window_size = 30, verbose = 0)
  mstomp_par_test1 <- mstomp_par(mp_toy_data$data[1:200, 1], window_size = 30, verbose = 0)
  # MSTOMP Multi
  mstomp_test <- mstomp(mp_toy_data$data[1:200, ], window_size = 30, verbose = 0)
  mstomp_test_must <- mstomp(mp_toy_data$data[1:200, ], window_size = 30, must_dim = c(1, 2), verbose = 0)
  mstomp_test_exc <- mstomp(mp_toy_data$data[1:200, ], window_size = 30, exc_dim = c(1, 2), verbose = 0)
  mstomp_par_test <- mstomp_par(mp_toy_data$data[1:200, ], window_size = 30, verbose = 0)
  mstomp_par_test_must <- mstomp_par(mp_toy_data$data[1:200, ], window_size = 30, must_dim = c(1, 2), verbose = 0)
  mstomp_par_test_exc <- mstomp_par(mp_toy_data$data[1:200, ], window_size = 30, exc_dim = c(1, 2), verbose = 0)

  if (skip_on_travis()) {
    test_that("Result hashes", {
      expect_known_hash(stamp_test, "1016c61c9f")
      expect_known_hash(stamp_join_test, "585be5fedc")
      expect_known_hash(stomp_test, "bd71cbd417")
      expect_known_hash(stomp_join_test, "5521be09db")
      expect_known_hash(mstomp_test, "fa0c150b92")
      expect_known_hash(mstomp_test1, "9c2bf3197d")
      expect_known_hash(mstomp_test_must, "13cefe2517")
      expect_known_hash(mstomp_test_exc, "b872f44cc5")
    })
  }

  # stamp_test and stamp_par_test
  test_that("Stamp equals to Stamp_par", {
    expect_equal(stamp_test, stamp_par_test)
  })

  # stamp_join_test and stamp_par_join_test
  test_that("Stamp Join equals to Stamp_par Join", {
    expect_equal(stamp_join_test, stamp_par_join_test)
  })

  # stamp_test and stomp_test
  test_that("Stamp equals to Stomp_par", {
    expect_equal(stamp_test, stomp_test)
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
