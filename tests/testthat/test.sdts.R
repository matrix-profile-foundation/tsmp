context("Testing SDTS functions")
library(tsmp)

w <- c(110, 220)
subs <- 11000:20000
tr_data <- test_data$train$data[subs]
tr_label <- test_data$train$label[subs]
te_data <- test_data$test$data[subs]
te_label <- test_data$test$label[subs]
model <- sdts.train(tr_data, tr_label, w, verbose = 0)
predict <- sdts.predict(model, te_data, round(mean(w)))
pred.score <- sdts.f.score(te_label, predict, 1)

test_that("SDTS Train", {
  expect_equal(round(model$score, 3) , 0.353)
  expect_equal(round(model$score.hist, 3), 0.353)
  expect_equal(round(sum(model$pattern[[1]])/sd(model$pattern[[1]]), 3), -4847.100)
  expect_equal(round(model$thold, 3), 8.726)
})

test_that("SDTS Predict", {
  expect_known_hash(predict, "1c340d263b")
})

test_that("SDTS F Score", {
  expect_equal(pred.score$f.score, 0.125)
  expect_equal(round(pred.score$precision, 4), 0.0667)
  expect_equal(pred.score$recall, 1)
})
