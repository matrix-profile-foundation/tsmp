context("Testing if Stamps and Stomps algorithms agree")
library(tsmp)

if (skip_on_cran()) {
  test_that("Errors", {
    # big window size
    expect_error(mstomp(toy_data$data[1:200, ], window.size = 500), regexp = "too short relative")
    expect_error(mstomp.par(toy_data$data[1:200, ], window.size = 500), regexp = "too short relative")
    expect_error(stomp(toy_data$data[1:200, ], window.size = 500), regexp = "too short relative")
    expect_error(stomp.par(toy_data$data[1:200, ], window.size = 500), regexp = "too short relative")
    expect_error(stamp(toy_data$data[1:200, ], window.size = 500), regexp = "too short relative")
    expect_error(stamp.par(toy_data$data[1:200, ], window.size = 500), regexp = "too short relative")

    # intersect
    expect_error(mstomp(toy_data$data[1:200, ], window.size = 30, must.dim = c(1, 2), exc.dim = c(2, 3)), regexp = "presented in both")
    expect_error(mstomp.par(toy_data$data[1:200, ], window.size = 30, must.dim = c(1, 2), exc.dim = c(2, 3)), regexp = "presented in both")
    # too many must.dim
    expect_error(mstomp(toy_data$data[1:200, ], window.size = 30, must.dim = c(1, 2, 3, 4)), regexp = "must have dimension must be less")
    expect_error(mstomp.par(toy_data$data[1:200, ], window.size = 30, must.dim = c(1, 2, 3, 4)), regexp = "must have dimension must be less")
    # too many exc.dim
    expect_error(mstomp(toy_data$data[1:200, ], window.size = 30, exc.dim = c(1, 2, 3, 4)), regexp = "exclusion dimension must be less")
    expect_error(mstomp.par(toy_data$data[1:200, ], window.size = 30, exc.dim = c(1, 2, 3, 4)), regexp = "exclusion dimension must be less")

    # small window size
    expect_error(stamp(toy_data$data[1:200, ], window.size = 2), regexp = "Subsequence length")
    expect_error(stamp.par(toy_data$data[1:200, ], window.size = 2), regexp = "Subsequence length")
    expect_error(mstomp(toy_data$data[1:200, ], window.size = 2), regexp = "Subsequence length")
    expect_error(mstomp.par(toy_data$data[1:200, ], window.size = 2), regexp = "Subsequence length")
    expect_error(stomp(toy_data$data[1:200, ], window.size = 2), regexp = "Subsequence length")
    expect_error(stomp.par(toy_data$data[1:200, ], window.size = 2), regexp = "Subsequence length")

    # unknown data type
    expect_error(stamp(table(rpois(100, 5)), window.size = 30), regexp = "Unknown type")
    expect_error(stamp.par(table(rpois(100, 5)), window.size = 30), regexp = "Unknown type")
    expect_error(mstomp(table(rpois(100, 5)), window.size = 30), regexp = "Unknown type")
    expect_error(mstomp.par(table(rpois(100, 5)), window.size = 30), regexp = "Unknown type")
    expect_error(stomp(table(rpois(100, 5)), window.size = 30), regexp = "Unknown type")
    expect_error(stomp.par(table(rpois(100, 5)), window.size = 30), regexp = "Unknown type")
  })

  test_that("Finish", {
    expect_message(stamp(toy_data$data[1:200, 1], window.size = 30), regex = "Finished")
    expect_message(stamp.par(toy_data$data[1:200, 1], window.size = 30), regex = "Finished")
    expect_message(mstomp(toy_data$data[1:200, 1], window.size = 30), regex = "Finished")
    expect_message(mstomp.par(toy_data$data[1:200, 1], window.size = 30), regex = "Finished")
    expect_message(stomp(toy_data$data[1:200, 1], window.size = 30), regex = "Finished")
    expect_message(stomp.par(toy_data$data[1:200, 1], window.size = 30), regex = "Finished")
  })

  # STAMP
  stamp.test <- stamp(toy_data$data[1:200, 1], window.size = 30, verbose = 0)
  stamp.join.test <- stamp(toy_data$data[1:200, 1], toy_data$data[1:100, 2], window.size = 30, verbose = 0)
  stamp.par.test <- stamp.par(toy_data$data[1:200, 1], window.size = 30, verbose = 0)
  stamp.par.join.test <- stamp.par(toy_data$data[1:200, 1], toy_data$data[1:100, 2], window.size = 30, verbose = 0)

  # STOMP
  stomp.test <- stomp(toy_data$data[1:200, 1], window.size = 30, verbose = 0)
  # stomp.join.test <- stomp(toy_data$data[1:200, 1], toy_data$data[1:100, 2], window.size = 30, verbose = 0)
  stomp.par.test <- stomp.par(toy_data$data[1:200, 1], window.size = 30, verbose = 0)
  # stomp.par.join.test <- stomp.par(toy_data$data[1:200, 1], toy_data$data[1:100, 2], window.size = 30, verbose = 0)

  # MSTOMP Uni
  mstomp.test1 <- mstomp(toy_data$data[1:200, 1], window.size = 30, verbose = 0)
  mstomp.par.test1 <- mstomp.par(toy_data$data[1:200, 1], window.size = 30, verbose = 0)
  # MSTOMP Multi
  mstomp.test <- mstomp(toy_data$data[1:200, ], window.size = 30, verbose = 0)
  mstomp.test.must <- mstomp(toy_data$data[1:200, ], window.size = 30, must.dim = c(1, 2), verbose = 0)
  mstomp.test.exc <- mstomp(toy_data$data[1:200, ], window.size = 30, exc.dim = c(1, 2), verbose = 0)
  mstomp.par.test <- mstomp.par(toy_data$data[1:200, ], window.size = 30, verbose = 0)
  mstomp.par.test.must <- mstomp.par(toy_data$data[1:200, ], window.size = 30, must.dim = c(1, 2), verbose = 0)
  mstomp.par.test.exc <- mstomp.par(toy_data$data[1:200, ], window.size = 30, exc.dim = c(1, 2), verbose = 0)

  if (skip_on_travis()) {
    test_that("Result hashes", {
      expect_known_hash(stamp.test, "1016c61c9f")
      expect_known_hash(stamp.join.test, "585be5fedc")
      expect_known_hash(stomp.test, "9c2bf3197d")
      # expect_known_hash(stomp.join.test, "585be5fedc")
      expect_known_hash(mstomp.test, "fa0c150b92")
      expect_known_hash(mstomp.test1, "9c2bf3197d")
      expect_known_hash(mstomp.test.must, "13cefe2517")
      expect_known_hash(mstomp.test.exc, "b872f44cc5")
    })
  }

  # stamp.test -> stamp.par.test
  # stamp.join.test -> stamp.par.join.test
  # stamp.test -> stomp.test
  # stamp.test -> stomp.par.test
  # stamp.test -> mstomp.test1
  # stamp.test -> mstomp.par.test1

  # mstomp.test -> mstomp.par.test
  # mstomp.test.must -> mstomp.par.test.must
  # mstomp.test.exc -> mstomp.par.test.exc

  # stamp.test -> stamp.par.test
  test_that("Stamp equals to Stamp.par", {
    expect_equal(stamp.test, stamp.par.test)
  })

  # stamp.join.test -> stamp.par.join.test
  test_that("Stamp Join equals to Stamp.par Join", {
    expect_equal(stamp.join.test, stamp.par.join.test)
  })

  # stamp.test -> stomp.test
  test_that("Stamp equals to Stomp.par", {
    expect_equal(stamp.test, stomp.test)
  })

  # test_that("Stomp equals to Stomp Join", {
  #   expect_equal(stomp.test, stomp.join.test)
  # })

  # test_that("Stomp Join equals to Stomp.par Join", {
  #   expect_equal(stomp.join.test, stomp.par.join.test)
  # })

  # stamp.test -> stomp.par.test
  test_that("Stamp equals to Stomp.par", {
    expect_equal(stamp.test, stomp.par.test)
  })

  # stamp.test -> mstomp.test1
  test_that("Stamp equals to mStomp", {
    expect_equal(stamp.test, mstomp.test1)
  })

  # stamp.test -> mstomp.par.test1
  test_that("Stamp equals to Stomp", {
    expect_equal(stamp.test, mstomp.par.test1)
  })


  # mstomp.test -> mstomp.par.test
  test_that("mStomp equals to mStomp.par", {
    expect_equal(mstomp.test, mstomp.par.test)
  })

  # mstomp.test.must -> mstomp.par.test.must
  test_that("mStomp must equals to mStomp.par must", {
    expect_equal(mstomp.test.must, mstomp.par.test.must)
  })

  # mstomp.test.exc -> mstomp.par.test.exc
  test_that("mStomp exc equals to mStomp.par exc", {
    expect_equal(mstomp.test.exc, mstomp.par.test.exc)
  })
}
