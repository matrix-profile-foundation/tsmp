test_that("Contrast Profile", {
  obj <- helper_contrast()
  cp <- contrast(obj$data1, obj$data2, obj$w, progress = FALSE)
  expect_snapshot_value(cp, style = "deparse")
})

test_that("Contrast Profile Parallel", {
  obj <- helper_contrast()
  cp <- contrast(obj$data1, obj$data2, obj$w, n_workers = 2L, progress = FALSE)
  expect_snapshot_value(cp, style = "deparse")
})

test_that("Results are equal", {
  obj <- helper_contrast()
  cp <- contrast(obj$data1, obj$data2, obj$w, n_workers = 1L, progress = FALSE)
  cp_par <- contrast(obj$data1, obj$data2, obj$w, n_workers = 2L, progress = FALSE)
  expect_equal(cp, cp_par)
})
