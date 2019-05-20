if (skip_on_cran()) {
  context("Testing Subset functions")
  library(tsmp)

  data <- mp_fluss_data$tilt_abp$data[20000:30000]
  w <- mp_fluss_data$tilt_abp$window
  nseg <- 1
  offset <- 4000
  test_mp <- tsmp(data, window_size = w, n_workers = 6)

  cac <- fluss_cac(test_mp)
  segments <- fluss_extract(cac, nseg)
  chain <- find_chains(test_mp)
  motif <- find_motif(test_mp)
  discord <- find_discord(test_mp)

  test_that("Sub Motif", {
    s_motif <- motif[1000:3000]
    expect_equal(s_motif$motif$motif_idx[[1]], c(366, 1147))
    expect_equal(s_motif$motif$motif_neighbor[[1]], c(1348, 560, 175, 1552, 951))
  })

  test_that("Sub Discord", {
    s_discord <- discord[1000:9000]
    expect_equal(s_discord$discord$discord_idx[[1]], 3891)
    expect_equal(s_discord$discord$discord_neighbor[[1]], 633)
  })

  test_that("Head Chain", {
    h_chain <- head(chain, 8000)
    expect_equal(sum(h_chain$chain$best), 39574)
    expect_equal(length(h_chain$chain$best), 6)
  })

  test_that("Tail Chain", {
    t_chain <- tail(chain, 4000)
    expect_equal(sum(t_chain$chain$best), 14625)
    expect_equal(length(t_chain$chain$best), 9)
  })

  test_that("Corrected Arc Count", {
    expect_equal(round(mean(cac$cac), 3), 0.377)
    expect_equal(round(sd(cac$cac), 3), 0.341)
    expect_equal(round(min(cac$cac), 3), 0)
    expect_equal(max(cac$cac), 1)
  })

  test_that("Head Arc Count", {
    h_cac <- head(cac, offset)
    expect_equal(round(mean(h_cac$cac), 4), 0.7198)
    expect_equal(round(sd(h_cac$cac), 3), 0.316)
    expect_equal(round(min(h_cac$cac), 3), 0.209)
    expect_equal(max(h_cac$cac), 1)
  })

  test_that("Tail Arc Count", {
    t_cac <- tail(cac, offset)
    expect_equal(round(mean(t_cac$cac), 4), 0.7654)
    expect_equal(round(sd(t_cac$cac), 3), 0.273)
    expect_equal(round(min(t_cac$cac), 3), 0.255)
    expect_equal(max(t_cac$cac), 1)
  })

  test_that("Segments found", {
    expect_equal(segments$fluss, 4902)
  })

  test_that("Head Segments found", {
    h_segment <- head(segments, offset)
    expect_equal(h_segment$fluss, 1228)
  })

  test_that("Tail Segments found", {
    t_segment <- tail(segments, offset)
    expect_equal(t_segment$fluss, 1643)
  })
}
