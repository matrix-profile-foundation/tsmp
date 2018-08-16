context("Testing if Stamps and Stomps algorithms agree")
library(tsmp)

if (skip_on_cran()) {
  test_that("Errors", {
    # big window size
    expect_error(mstomp(toy_data$data[1:200, ], window.size = 500, verbose = 0), regexp = "too short relative")
    expect_error(mstomp.par(toy_data$data[1:200, ], window.size = 500, verbose = 0), regexp = "too short relative")
    expect_error(stamp(toy_data$data[1:200, ], window.size = 500, verbose = 0), regexp = "too short relative")
    expect_error(stamp.par(toy_data$data[1:200, ], window.size = 500, verbose = 0), regexp = "too short relative")

    # intersect
    expect_error(mstomp(toy_data$data[1:200, ], window.size = 30, must.dim = c(1, 2), exc.dim = c(2, 3), verbose = 0), regexp = "presented in both")
    expect_error(mstomp.par(toy_data$data[1:200, ], window.size = 30, must.dim = c(1, 2), exc.dim = c(2, 3), verbose = 0), regexp = "presented in both")
    # too many must.dim
    expect_error(mstomp(toy_data$data[1:200, ], window.size = 30, must.dim = c(1, 2, 3, 4), verbose = 0), regexp = "must have dimension must be less")
    expect_error(mstomp.par(toy_data$data[1:200, ], window.size = 30, must.dim = c(1, 2, 3, 4), verbose = 0), regexp = "must have dimension must be less")
    # too many exc.dim
    expect_error(mstomp(toy_data$data[1:200, ], window.size = 30, exc.dim = c(1, 2, 3, 4), verbose = 0), regexp = "exclusion dimension must be less")
    expect_error(mstomp.par(toy_data$data[1:200, ], window.size = 30, exc.dim = c(1, 2, 3, 4), verbose = 0), regexp = "exclusion dimension must be less")

    # small window size
    expect_error(stamp(toy_data$data[1:200, ], window.size = 2, verbose = 0), regexp = "Subsequence length")
    expect_error(stamp.par(toy_data$data[1:200, ], window.size = 2, verbose = 0), regexp = "Subsequence length")
    expect_error(mstomp(toy_data$data[1:200, ], window.size = 2, verbose = 0), regexp = "Subsequence length")
    expect_error(mstomp.par(toy_data$data[1:200, ], window.size = 2, verbose = 0), regexp = "Subsequence length")

    # unknown data type
    expect_error(stamp(table(rpois(100, 5)), window.size = 30, verbose = 0), regexp = "Unknown type")
    expect_error(stamp.par(table(rpois(100, 5)), window.size = 30, verbose = 0), regexp = "Unknown type")
    expect_error(mstomp(table(rpois(100, 5)), window.size = 30, verbose = 0), regexp = "Unknown type")
    expect_error(mstomp.par(table(rpois(100, 5)), window.size = 30, verbose = 0), regexp = "Unknown type")
  })

  stamp.test <- stamp(toy_data$data[1:200, 1], window.size = 30, verbose = 0)
  Sys.sleep(0.5)
  stamp.join.test <- stamp(toy_data$data[1:200, 1], toy_data$data[1:100, 2], window.size = 30, verbose = 0)
  Sys.sleep(0.5)
  stamp.par.test <- stamp.par(toy_data$data[1:200, 1], window.size = 30, verbose = 0)
  Sys.sleep(0.5)
  stamp.par.join.test <- stamp.par(toy_data$data[1:200, 1], toy_data$data[1:100, 2], window.size = 30, verbose = 0)
  Sys.sleep(0.5)
  stomp.test <- mstomp(toy_data$data[1:200, 1], window.size = 30, verbose = 0)
  Sys.sleep(0.5)
  stomp.par.test <- mstomp.par(toy_data$data[1:200, 1], window.size = 30, verbose = 0)
  Sys.sleep(0.5)
  mstomp.test <- mstomp(toy_data$data[1:200, ], window.size = 30, verbose = 0)
  Sys.sleep(0.5)
  mstomp.test.must <- mstomp(toy_data$data[1:200, ], window.size = 30, must.dim = c(1, 2), verbose = 0)
  Sys.sleep(0.5)
  mstomp.test.exc <- mstomp(toy_data$data[1:200, ], window.size = 30, exc.dim = c(1, 2), verbose = 0)
  Sys.sleep(0.5)
  mstomp.par.test <- mstomp.par(toy_data$data[1:200, ], window.size = 30, verbose = 0)
  Sys.sleep(0.5)
  mstomp.par.test.must <- mstomp.par(toy_data$data[1:200, ], window.size = 30, must.dim = c(1, 2), verbose = 0)
  Sys.sleep(0.5)
  mstomp.par.test.exc <- mstomp.par(toy_data$data[1:200, ], window.size = 30, exc.dim = c(1, 2), verbose = 0)

  test_that("Result hashes", {
    expect_known_hash(stamp.test, "1016c61c9f")
    expect_known_hash(stamp.join.test, "c644733f25")
    expect_known_hash(mstomp.test, "fa0c150b92")
    expect_known_hash(mstomp.test.must, "13cefe2517")
    expect_known_hash(mstomp.test.exc, "b872f44cc5")
  })

  test_that("Stamp equals to Stamp.par", {
    expect_equal(stamp.test, stamp.par.test)
  })

  test_that("Stamp Join equals to Stamp.par Join", {
    expect_equal(stamp.join.test, stamp.par.join.test)
  })

  test_that("Stomp equals to Stomp.par", {
    expect_equal(stomp.test, stomp.par.test)
  })

  test_that("Stamp equals to Stomp", {
    expect_equal(stamp.test, stomp.test)
  })

  test_that("Stamp equals to Stomp.par", {
    expect_equal(stamp.test, stomp.par.test)
  })

  test_that("Stomp equals to Stamp.par", {
    expect_equal(stomp.test, stamp.par.test)
  })

  test_that("mStomp equals to mStomp.par", {
    expect_equal(mstomp.test, mstomp.par.test)
  })

  test_that("mStomp must equals to mStomp.par must", {
    expect_equal(mstomp.test.must, mstomp.par.test.must)
  })

  test_that("mStomp exc equals to mStomp.par exc", {
    expect_equal(mstomp.test.exc, mstomp.par.test.exc)
  })
}
