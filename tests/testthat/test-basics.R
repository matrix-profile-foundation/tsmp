if (skip_on_cran()) {
  context("Testing if basic functions are ok")
  library(tsmp)
  w <- 30
  ref_data <- mp_toy_data$data[, 1]
  query_data <- mp_toy_data$data[, 1]
  d_size <- length(ref_data)
  q_size <- length(query_data)

  test_that("Errors", {
    # big window size
    expect_error(fast_movsd(mp_toy_data$data[, 1], 1), "must be at least 2")
    expect_error(fast_movsd(mp_toy_data$data[1:100, 1], 500), "is too large")
    expect_error(expect_message(beep(audio::close.audioInstance(99)), "Failed"))
    expect_error(tsmp:::diff2(data.frame(1:10), as.matrix(10:1)), "matrices")
    expect_error(tsmp:::diff2(as.matrix(1:10), matrix(10:1, ncol = 2)), "columns")
  })

  pre <- mass_pre(ref_data, query_data, w)
  pre$query_mean <- pre$query_mean[1]
  pre$query_sd <- pre$query_sd[1]

  pre3 <- c(pre, list(data = ref_data, k = NULL))

  res <- do.call("mass_v2", (c(list(query_data[1:w]), pre)))
  res3 <- do.call("mass_v3", (c(list(query_data[1:w]), pre3)))
  pre3$k <- 1 # fix hash later

  movsd <- fast_movsd(mp_toy_data$data[, 1], 30)
  movavg <- fast_movavg(mp_toy_data$data[, 1], 30)

  prew <- mass_pre_w(ref_data, query_data, w, c(rep(1, 15), rep(0.5, 15)))
  resw <- do.call("mass_weighted", (c(list(query_data[1:w]), prew)))



  test_that("Fast Moving SD is ok", {
    expect_known_hash(round(movsd, 3), "ffda40fd35")
  })

  test_that("Fast Moving Average is ok", {
    expect_known_hash(round(movavg, 3), "601febe569")
  })

  test_that("MASS Pre is ok", {
    expect_known_hash(lapply(pre, round, 3), "1ac5ece94f")
    expect_known_hash(lapply(pre3, round, 3), "3fbcbbbc6f")
    expect_known_hash(lapply(prew, round, 3), "7199ea592a")
  })

  test_that("MASS is ok", {
    expect_equal(sum(round(Re(res$distance_profile), 2)), 30737.17)
    expect_equal(sum(round(Re(res$last_product), 2)), 5965.13)
    expect_equal(sum(round(Re(res3$distance_profile), 2)), 30737.17)
    expect_equal(sum(round(Re(res3$last_product), 2)), 5965.13)
    expect_equal(sum(round(Re(resw$distance_profile), 2)), 21944.89)
    expect_equal(sum(round(Re(resw$last_product), 2)), -1076.72)
  })

  pred <- rep(c(rep(1, 220), rep(0, 7830)), 27)[1:214791]
  fs1 <- sdts_score(pred, mp_test_data$test$label, 1)
  fs2 <- sdts_score(pred, mp_test_data$train$label, 1)

  test_that("F-Score is ok", {
    expect_known_hash(lapply(fs1, round, 3), "8241d8411e")
    expect_known_hash(lapply(fs2, round, 3), "ad8b052f2c")
  })
}
