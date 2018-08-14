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
  expect_known_hash(round(movsd, 3), "ffda40fd35")
})

test_that("Fast Moving Average is ok", {
  expect_known_hash(round(movavg, 3), "601febe569")
})

test_that("MASS Pre is ok", {
  expect_known_hash(lapply(pre, round, 3), "247c3e49cd")
})

test_that("MASS is ok", {
  expect_equal(sum(round(Re(res$distance.profile), 2)), 30737.17)
  expect_equal(sum(round(Re(res$last.product), 2)), 5965.13)
})

pred <- rep(c(rep(1, 220), rep(0, 7830)), 27)[1:214791]
fs1 <- sdts.f.score(test_data$test$label, pred, 1)
fs2 <- sdts.f.score(test_data$train$label, pred, 1)

test_that("F-Score is ok", {
  expect_known_hash(lapply(fs1, round, 3), "d016bdcb98")
  expect_known_hash(lapply(fs2, round, 3), "7b1492f832")
})
