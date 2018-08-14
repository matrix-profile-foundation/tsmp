context("Testing if Stamps and Stomps algorithms agree")
library(tsmp)

if (skip_on_cran()) {
  stamp.test <- stamp(toy_data$data[1:200, 1], window.size = 30, verbose = 0)
  Sys.sleep(1)
  stamp.par.test <- stamp.par(toy_data$data[1:200, 1], window.size = 30, verbose = 0)
  Sys.sleep(1)
  stomp.test <- mstomp(toy_data$data[1:200, 1], window.size = 30, verbose = 0)
  Sys.sleep(1)
  stomp.par.test <- mstomp.par(toy_data$data[1:200, 1], window.size = 30, verbose = 0)
  Sys.sleep(1)
  mstomp.test <- mstomp(toy_data$data[1:200, ], window.size = 30, verbose = 0)
  Sys.sleep(1)
  mstomp.par.test <- mstomp.par(toy_data$data[1:200, ], window.size = 30, verbose = 0)

  test_that("Stamp equals to Stamp.par", {
    expect_equal(stamp.test, stamp.par.test)
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
}
