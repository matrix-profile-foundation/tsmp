test_that("as_ts_data handles NULL and empty inputs", {
  expect_null(as_ts_data(NULL))
  expect_null(convert_data(NULL))
  expect_equal(as_ts_data(list()), numeric(0))
})

test_that("as_ts_data standardizes 1D univariate inputs to atomic vectors", {
  # Atomic numeric and integer vectors
  v1 <- 1:10
  expect_identical(as_ts_data(v1), as.numeric(1:10))
  expect_null(dim(as_ts_data(v1)))

  # 1D Column matrix (N x 1)
  m_col <- matrix(1:10, ncol = 1)
  res_col <- as_ts_data(m_col)
  expect_identical(res_col, as.numeric(1:10))
  expect_null(dim(res_col))

  # 1D Row matrix (1 x N)
  m_row <- matrix(1:10, nrow = 1)
  res_row <- as_ts_data(m_row)
  expect_identical(res_row, as.numeric(1:10))
  expect_null(dim(res_row))

  # 1-column data.frame
  df_1 <- data.frame(val = 1:10)
  res_df1 <- as_ts_data(df_1)
  expect_identical(res_df1, as.numeric(1:10))
  expect_null(dim(res_df1))

  # 1-element list
  lst_1 <- list(1:10)
  res_lst1 <- as_ts_data(lst_1)
  expect_identical(res_lst1, as.numeric(1:10))
  expect_null(dim(res_lst1))
})

test_that("as_ts_data standardizes multivariate 2D inputs to N x D matrices", {
  # Vertical 2D matrix (N >= D)
  m_vert <- matrix(1:20, nrow = 10, ncol = 2)
  res_vert <- as_ts_data(m_vert)
  expect_true(is.matrix(res_vert))
  expect_equal(dim(res_vert), c(10, 2))
  expect_equal(storage.mode(res_vert), "double")

  # Horizontal 2D matrix (rows < cols) -> transposed to vertical (N x D)
  m_horiz <- matrix(1:20, nrow = 2, ncol = 10)
  res_horiz <- as_ts_data(m_horiz)
  expect_true(is.matrix(res_horiz))
  expect_equal(dim(res_horiz), c(10, 2))
  expect_equal(storage.mode(res_horiz), "double")

  # Multi-column data frame
  df_multi <- data.frame(a = 1:10, b = 11:20)
  res_df_multi <- as_ts_data(df_multi)
  expect_true(is.matrix(res_df_multi))
  expect_equal(dim(res_df_multi), c(10, 2))
  expect_equal(storage.mode(res_df_multi), "double")
})

test_that("as_ts_data standardizes lists of homogeneous series", {
  # List of two 1D vectors of equal length
  l_vecs <- list(1:10, 11:20)
  res_l_vecs <- as_ts_data(l_vecs)
  expect_true(is.matrix(res_l_vecs))
  expect_equal(dim(res_l_vecs), c(10, 2))
  expect_equal(res_l_vecs[, 1], as.numeric(1:10))
  expect_equal(res_l_vecs[, 2], as.numeric(11:20))

  # List of 1D matrices
  l_mats <- list(matrix(1:10, ncol = 1), matrix(11:20, nrow = 1))
  res_l_mats <- as_ts_data(l_mats)
  expect_true(is.matrix(res_l_mats))
  expect_equal(dim(res_l_mats), c(10, 2))
})

test_that("as_ts_data raises appropriate errors for invalid structures", {
  # 3D arrays
  arr3d <- array(1:24, dim = c(2, 3, 4))
  expect_error(as_ts_data(arr3d), "cannot have more than 2 dimensions")

  # List of series with different lengths
  l_diff <- list(1:10, 1:5)
  expect_error(as_ts_data(l_diff), "must have the same length")

  # List containing multivariate matrices
  l_multi_mat <- list(1:10, matrix(1:20, nrow = 10, ncol = 2))
  expect_error(as_ts_data(l_multi_mat), "must be univariate")
})
