if (!testthat:::on_cran()) {
  context("Testing Annotation functions")
  library(tsmp)

  data <- as.matrix(mp_test_data$train$data[1:1000])
  window <- 50
  mp <- tsmp(data, window_size = window, verbose = 0)
  stop_loc <- 150
  prof_size <- nrow(data) - window + 1
  comp <- hard <- motion <- stopw <- zero <- NULL
  avcomp <- av_complexity(mp, apply = TRUE)

  test_that("Silent", {
    expect_silent(comp <<- av_complexity(mp))
    expect_silent(hard <<- av_hardlimit_artifact(mp))
    expect_silent(motion <<- av_motion_artifact(mp))
    expect_silent(stopw <<- av_stop_word(mp, stop_word_loc = stop_loc))
    expect_silent(zero <<- av_zerocrossing(mp))
  })

  test_that("Apply", {
    expect_equal(av_apply(comp), avcomp)
    expect_error(av_apply(avcomp), "already")
    class(avcomp) <- "MatrixProfile"
    expect_error(av_apply(avcomp), "class `AnnotationVector`")
    class(avcomp) <- "test"
    expect_error(av_apply(avcomp), "class `MatrixProfile`")
  })

  test_that("Result dim", {
    expect_equal(dim(comp$av), c(prof_size, 1))
    expect_equal(dim(hard$av), c(prof_size, 1))
    expect_equal(dim(motion$av), c(prof_size, 1))
    expect_equal(dim(stopw$av), c(prof_size, 1))
    expect_equal(dim(zero$av), c(prof_size, 1))
  })
  test_that("Result values", {
    expect_equal(round(sum(comp$av) / sd(comp$av), 2), 1689.92)
    expect_equal(round(sum(hard$av) / sd(hard$av), 2), 3568.52)
    expect_equal(round(sum(motion$av) / sd(motion$av), 1), 1015.7)
    expect_equal(round(sum(stopw$av) / sd(stopw$av), 2), 1336.86)
    expect_equal(round(sum(zero$av) / sd(zero$av), 2), 666.75)
  })
}
