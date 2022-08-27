test_that("Basic distance profile", {
  withr::local_options(digits = 22)
  obj <- helper_dist_profile()

  nn <- dist_profile(data = obj$ref_data, window_size = obj$w)
  dp <- sqrt(nn$distance_profile)

  nn %>%
    expect_type("list") %>%
    expect_length(3) %>%
    expect_named(c("distance_profile", "last_product", "params"))

  expect_equal(nn$distance_profile, dp * dp)
  expect_snapshot(nn)
})

test_that("Non-normalized distance profile", {
  withr::local_options(digits = 22)
  obj <- helper_dist_profile()

  nn <- dist_profile(data = obj$ref_data, window_size = obj$w, type = "non_normalized")
  dp <- sqrt(nn$distance_profile)

  nn %>%
    expect_type("list") %>%
    expect_length(3) %>%
    expect_named(c("distance_profile", "last_product", "params"))

  expect_equal(nn$distance_profile, dp * dp)
  expect_snapshot(nn)
})

test_that("Absolute distance profile", {
  withr::local_options(digits = 22)
  obj <- helper_dist_profile()

  nn <- dist_profile(data = obj$ref_data, window_size = obj$w, type = "absolute")
  dp <- sqrt(nn$distance_profile)

  nn %>%
    expect_type("list") %>%
    expect_length(3) %>%
    expect_named(c("distance_profile", "last_product", "params"))

  expect_equal(nn$distance_profile, dp * dp)
  expect_snapshot(nn)
})

test_that("Weighted distance profile", {
  withr::local_options(digits = 22)
  obj <- helper_dist_profile()

  weights <- abs(sin(seq(0, 3, length.out = obj$w) * 3))
  nn <- dist_profile(
    data = obj$ref_data, window_size = obj$w, type = "weighted",
    weights = weights
  )
  dp <- sqrt(nn$distance_profile)

  nn %>%
    expect_type("list") %>%
    expect_length(3) %>%
    expect_named(c("distance_profile", "last_product", "params"))

  expect_equal(nn$distance_profile, dp * dp)
  expect_snapshot(nn)
})

test_that("Basic distance profile AB", {
  withr::local_options(digits = 22)
  obj <- helper_dist_profile()

  nn <- dist_profile(data = obj$ref_data, window_size = obj$w, query = obj$query_data)
  dp <- sqrt(nn$distance_profile)

  nn %>%
    expect_type("list") %>%
    expect_length(3) %>%
    expect_named(c("distance_profile", "last_product", "params"))

  expect_equal(nn$distance_profile, dp * dp)
  expect_snapshot(nn)
})

test_that("Non-normalized distance profile AB", {
  withr::local_options(digits = 22)
  obj <- helper_dist_profile()

  nn <- dist_profile(data = obj$ref_data, window_size = obj$w, query = obj$query_data, type = "non_normalized")
  dp <- sqrt(nn$distance_profile)

  nn %>%
    expect_type("list") %>%
    expect_length(3) %>%
    expect_named(c("distance_profile", "last_product", "params"))

  expect_equal(nn$distance_profile, dp * dp)
  expect_snapshot(nn)
})

test_that("Absolute distance profile AB", {
  withr::local_options(digits = 22)
  obj <- helper_dist_profile()

  nn <- dist_profile(data = obj$ref_data, window_size = obj$w, query = obj$query_data, type = "absolute")
  dp <- sqrt(nn$distance_profile)

  nn %>%
    expect_type("list") %>%
    expect_length(3) %>%
    expect_named(c("distance_profile", "last_product", "params"))

  expect_equal(nn$distance_profile, dp * dp)
  expect_snapshot(nn)
})

test_that("Weighted distance profile AB", {
  withr::local_options(digits = 22)
  obj <- helper_dist_profile()

  weights <- abs(sin(seq(0, 3, length.out = obj$w) * 3))
  nn <- dist_profile(
    data = obj$ref_data, window_size = obj$w, type = "weighted", query = obj$query_data,
    weights = weights
  )
  dp <- sqrt(nn$distance_profile)

  nn %>%
    expect_type("list") %>%
    expect_length(3) %>%
    expect_named(c("distance_profile", "last_product", "params"))

  expect_equal(nn$distance_profile, dp * dp)
  expect_snapshot(nn)
})

test_that("Iterative distance profile", {
  withr::local_options(digits = 22)
  obj <- helper_dist_profile()

  nn <- dist_profile(data = obj$ref_data, window_size = obj$w)
  nn2 <- dist_profile(data = obj$ref_data, window_size = obj$w, index = 10)
  nn3 <- dist_profile(data = obj$ref_data, window_size = obj$w, params = nn, index = 10)
  expect_equal(nn2, nn3)
})
