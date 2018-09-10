if (skip_on_cran()) {
  context("Testing Salient functions")
  library(tsmp)

  # fake data for now
  w <- 100
  len <- 5000
  data <- as.matrix(mp_test_data$train$data[1:len])
  mps <- (5000 - w + 1)
  mp <- (data - mean(data))[1:mps]
  pi <- floor((rev(mp) - min(mp)) / (max(mp) - min(mp)) * mps)
  mp <- as.matrix(mp)
  pi <- as.matrix(pi)
  res <- NULL

  dummy <- list()
  dummy$mp <- mp
  dummy$pi <- pi
  dummy$w <- w
  dummy$ez <- 0.5
  dummy$data <- list()
  dummy$data[[1]] <- data
  class(dummy) <- "MatrixProfile"

  test_that("Finish", {
    expect_silent(salient_subsequences(dummy, verbose = 0))

    if (skip_on_cran()) {
      expect_message(res <<- salient_subsequences(dummy, verbose = 2), regex = "Finished")
    } else {
      expect_message(res <<- salient_subsequences(dummy, verbose = 1), regex = "Finished")
    }
  })

  if (!is.null(res)) {
    test_that("Score", {
      expect_equal(round(sum(res$salient$indexes) / sd(res$salient$indexes), 4), 140.5346)
      expect_equal(round(sum(res$salient$idx_bit_size) / sd(res$salient$idx_bit_size), 2), 24166.32)
      expect_equal(tsmp:::get_bitsize(data > 0, 10), 2270)
      expect_equal(sum(tsmp:::discrete_norm(data, 3, max(data), min(data))), 36632)
      maxmin <- tsmp:::discrete_norm_pre(as.vector(data), w)
      expect_equal(round(maxmin$max, 4), 4.1591)
      expect_equal(round(maxmin$min, 4), -4.6241)
      expect_equal(tsmp:::get_sorted_idx(mp, 10), c(3761, 3995, 3996, 4001, 3997, 3835, 4002, 4007, 3832, 4005))
      expect_equal(round(sd(tsmp:::salient_mds(data, list(indexes = c(1200, 1400, 1600, 1800)), w)), 4), 4.3979)
      res$salient$idx_bit_size <- rev(res$salient$idx_bit_size)
      res$salient$idx_bit_size[45] <- 0
      truth <- res$salient$indexes
      truth[10:15] <- 0
      truth[20:30] <- 0
      truth <- truth + 10
      scr <- tsmp:::salient_score(truth, res$salient, w)
      expect_equal(round(scr$precision, 4), 0.6136)
      expect_equal(round(scr$recall, 4), 0.54)
      expect_equal(round(scr$fscore, 4), 0.5745)
    })
  }
}
