if (!testthat:::on_cran()) {
  context("Testing SDTS functions")
  library(tsmp)

  test_that("Errors", {
    # big window size
    expect_error(sdts_train(mp_test_data$train$data[1:100], mp_test_data$train$label[1:110],
      window_size = 5000
    ), "Time series is too short")

    # small window size
    expect_error(sdts_train(mp_test_data$train$data[1:100], mp_test_data$train$label[1:100],
      window_size = 2
    ), "window_size")

    # unknown data type
    expect_error(sdts_train(table(mp_test_data$train$data[1:100]), mp_test_data$train$label[1:100],
      window_size = 110
    ), "Unknown type")
  })

  w <- c(110, 220, 330)
  subs <- 20000:60000
  tr_data <- as.data.frame(mp_test_data$train$data[subs])
  tr_label <- mp_test_data$train$label[subs]
  te_data <- mp_test_data$test$data[subs]
  te_label <- mp_test_data$test$label[subs]
  model <- sdts_train(tr_data, tr_label, w, verbose = 0)
  predict <- sdts_predict(model, te_data, round(mean(w)))
  pred_score <- sdts_score(predict, te_label, 1)

  test_that("SDTS Train", {
    expect_equal(round(model$score, 3), 0.889)
    expect_equal(round(model$score_hist, 3), c(0.667, 0.889))
    expect_equal(round(sum(model$pattern[[1]] + model$pattern[[2]]) / sd(model$pattern[[1]]), 3), -8289.256)
    expect_equal(round(model$thold, 3), c(9.125, 2.069))
  })

  test_that("SDTS Predict", {
    expect_known_hash(predict, "72ddd7b33b")
  })

  test_that("SDTS F Score", {
    expect_equal(pred_score$f_score, 0.8)
    expect_equal(round(pred_score$precision, 4), 0.8)
    expect_equal(pred_score$recall, 0.8)
  })
}
