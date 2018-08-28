context("Testing SDTS functions")
library(tsmp)

test_that("Errors", {
  # big window size
  expect_error(sdts.train(test_data$train$data[1:100], test_data$train$label[1:110], window.size = 5000), regexp = "Time series is too short")

  # small window size
  expect_error(sdts.train(test_data$train$data[1:100], test_data$train$label[1:100], window.size = 2), regexp = "Subsequence length must")

  # unknown data type
  expect_error(sdts.train(table(test_data$train$data[1:100]), test_data$train$label[1:100], window.size = 110), regexp = "Unknown type")
})

w <- c(110, 220, 330)
subs <- 20000:60000
tr.data <- as.data.frame(test_data$train$data[subs])
tr.label <- test_data$train$label[subs]
te.data <- test_data$test$data[subs]
te.label <- test_data$test$label[subs]
model <- sdts.train(tr.data, tr.label, w, verbose = 0)
predict <- sdts.predict(model, te.data, round(mean(w)))
pred.score <- sdts.f.score(te.label, predict, 1)

test_that("SDTS Train", {
  expect_equal(round(model$score, 3), 0.889)
  expect_equal(round(model$score.hist, 3), c(0.667, 0.889))
  expect_equal(round(sum(model$pattern[[1]] + model$pattern[[2]]) / sd(model$pattern[[1]]), 3), -8289.256)
  expect_equal(round(model$thold, 3), c(9.125, 2.069))
})

test_that("SDTS Predict", {
  expect_known_hash(predict, "72ddd7b33b")
})

test_that("SDTS F Score", {
  expect_equal(pred.score$f.score, 0.8)
  expect_equal(round(pred.score$precision, 4), 0.8)
  expect_equal(pred.score$recall, 0.8)
})
