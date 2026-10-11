test_that("tsmp handles NA data with mpx and stomp", {
  data("mp_toy_data", package = "tsmp", envir = environment())
  data <- mp_toy_data$data[1:200, 1]
  data_na <- data
  data_na[50:55] <- NA_real_

  # 1. mpx with NA
  expect_no_error({
    mp_mpx <- tsmp(data_na, window_size = 30, mode = "mpx", verbose = 0)
  })
  expect_s3_class(mp_mpx, "MatrixProfile")
  expect_true(any(is.infinite(mp_mpx$mp) | is.na(mp_mpx$mp)))

  # 2. stomp with NA
  expect_no_error({
    mp_stomp <- tsmp(data_na, window_size = 30, mode = "stomp", verbose = 0)
  })
  expect_s3_class(mp_stomp, "MatrixProfile")
  expect_true(any(is.infinite(mp_stomp$mp)))

  # 3. find_motif on MP with NA
  expect_no_error({
    motifs <- find_motif(mp_mpx, n_motifs = 1, n_neighbors = 2)
  })
  expect_true(!is.null(motifs$motif))

  # 4. find_discord on MP with NA (must not pick Inf/NA as discord)
  expect_no_error({
    discords <- find_discord(mp_mpx, n_discords = 1)
  })
  expect_true(!is.null(discords$discord))
  d_idx <- discords$discord$discord_idx[[1]]
  expect_true(is.finite(mp_mpx$mp[d_idx]))
})

test_that("dist_profile supports n_workers and compute supports multi-core", {
  data("mp_toy_data", package = "tsmp", envir = environment())
  ref_data <- mp_toy_data$data[1:200, 1]
  
  # dist_profile with n_workers
  dp1 <- dist_profile(ref_data, ref_data, window_size = 30, index = 10, n_workers = 1)
  dp2 <- dist_profile(ref_data, ref_data, window_size = 30, index = 10, n_workers = 2)
  expect_equal(dp1$distance_profile, dp2$distance_profile)

  # compute with n_jobs > 1
  expect_no_error({
    cmp <- compute(ref_data, windows = 30, n_jobs = 2L)
  })
  expect_s3_class(cmp, "MatrixProfile")

  # find_snippet with n_workers
  expect_no_error({
    snip <- find_snippet(ref_data, s_size = 40, n_snippets = 2, n_workers = 2L)
  })
  expect_s3_class(snip, "Snippet")

  # pmp with n_workers
  expect_no_error({
    pmp_res <- pmp(ref_data, window_sizes = c(30, 40), n_workers = 2L, verbose = 0)
  })
  expect_s3_class(pmp_res, "PMP")
})
