context("Testing Salient functions")
library(tsmp)

# fake data for now
w <- 100
len <- 5000
data <- as.matrix(mp_test_data$train$data[1:len])
mps <- (length(data) - w + 1)
mp <- (data - mean(data))[1:mps]
pi <- floor((rev(mp) - min(mp)) / (max(mp) - min(mp)) * mps)
mp <- as.matrix(mp)
pi <- as.matrix(pi)
res <- NULL

test_that("Errors", {
  # unknown data type
  expect_error(salient_subsequences(table(data), matrix_profile = mp, profile_index = pi, window_size = w), regexp = "Unknown type")
})

test_that("Finish", {
  expect_silent(salient_subsequences(as.data.frame(t(data)), matrix_profile = mp, profile_index = pi, window_size = w, verbose = 0))
  expect_silent(salient_subsequences(as.vector(data), matrix_profile = mp, profile_index = pi, window_size = w, verbose = 0))
  expect_silent(salient_subsequences(list(as.vector(data)), matrix_profile = mp, profile_index = pi, window_size = w, verbose = 0))
  if (skip_on_cran()) {
    expect_message(res <<- salient_subsequences(data, matrix_profile = mp, profile_index = pi, window_size = w, verbose = 2), regex = "Finished")
  } else {
    expect_message(res <<- salient_subsequences(data, matrix_profile = mp, profile_index = pi, window_size = w, verbose = 1), regex = "Finished")
  }
})

if (!is.null(res)) {
  test_that("Score", {
    expect_equal(round(sum(res$indexes) / sd(res$indexes), 4), 140.5346)
    expect_equal(round(sum(res$idx_bit_size) / sd(res$idx_bit_size), 2), 24166.32)
    expect_equal(tsmp:::get_bitsize(data > 0, 10), 2270)
    expect_equal(sum(tsmp:::discrete_norm(data, 3, max(data), min(data))), 36632)
    maxmin <- tsmp:::discrete_norm_pre(as.vector(data), w)
    expect_equal(round(maxmin$max, 4), 4.1591)
    expect_equal(round(maxmin$min, 4), -4.6241)
    expect_equal(tsmp:::get_sorted_idx(mp, 10), c(3761, 3995, 3996, 4001, 3997, 3835, 4002, 4007, 3832, 4005))
    expect_equal(round(sd(tsmp:::salient_mds(data, list(indexes = c(1200, 1400, 1600, 1800)), w)), 4), 4.3979)
    res$idx_bit_size <- rev(res$idx_bit_size)
    res$idx_bit_size[45] <- 0
    truth <- res$indexes
    truth[10:15] <- 0
    truth[20:30] <- 0
    truth <- truth + 10
    scr <- tsmp:::salient_score(truth, res, w)
    expect_equal(round(scr$precision, 4), 0.6136)
    expect_equal(round(scr$recall, 4), 0.54)
    expect_equal(round(scr$fscore, 4), 0.5745)
  })
}
