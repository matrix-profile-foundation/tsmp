context("Testing if basic functions are ok")
library(tsmp)
w <- 30
ref.data <- toy_data$data[,1]
query.data <- toy_data$data[,1]
d.size <- length(ref.data)
q.size <- length(query.data)

pre <- mass.pre(ref.data, d.size, query.data, q.size, w)
res <- mass(pre$data.fft, query.data[1:w], d.size, w, pre$data.mean, pre$data.sd, pre$query.mean[1], pre$query.sd[1])
movsd <- fast.movsd(toy_data$data[,1], 30)
movavg <- fast.movavg(toy_data$data[,1], 30)

test_that("Fast Moving SD is ok", {
  expect_known_hash(movsd, "2b855c60d0")
})

test_that("Fast Moving Average is ok", {
  expect_known_hash(movavg, "141e89804f")
})

test_that("MASS Pre is ok", {
  expect_known_hash(pre, "de3fa83ee1")
})

test_that("MASS is ok", {
  expect_known_hash(res, "262c7d837a")
})

pred <- rep(c(rep(1, 220), rep(0, 7830)), 27)[1:214791]
fs1 <- sdts.f.score(test_data$test$label, pred, 1)
fs2 <- sdts.f.score(test_data$train$label, pred, 1)

test_that("F-Score is ok", {
  expect_known_hash(fs1, "74f2c6d1e0")
  expect_known_hash(fs2, "0d298030ea")
})
