if (!testthat:::on_cran()) {
  context("Testing Subset functions")
  library(tsmp)

  data <- mp_fluss_data$tilt_abp$data[20000L:30000L]
  mpd_query <- mp_fluss_data$tilt_abp$data[25001L:25210L]
  data_new <- mp_fluss_data$tilt_abp$data[30001L:31000L]
  w <- mp_fluss_data$tilt_abp$window
  nseg <- 1
  offset <- 4000
  test_mp <- tsmp(data, window_size = w)
  mpd_vect <- mpdist(data, mpd_query, floor(w / 10), type = "vector")
  test_floss <- floss(test_mp, data_new, 5000)

  cac <- fluss_cac(test_mp)
  segments <- fluss_extract(cac, nseg)
  chain <- find_chains(test_mp)
  motif <- find_motif(test_mp)
  discord <- find_discord(test_mp)
  salient <- salient_subsequences(test_mp, n_bits = c(4, 6, 8), verbose = 0)
  a_vector <- av_complexity(test_mp, apply = TRUE)

  test_that("Sub FLOSS", {
    s_test_floss <- test_floss[50L:4000L]
    expect_equal(round(sum(s_test_floss$cac) / sd(s_test_floss$cac), 1), 9316.1)
  })

  test_that("Sub MPDist", {
    s_mpd_vect <- mpd_vect[1000L:3000L]

    expect_equal(round(sum(s_mpd_vect$mpdist) / sd(s_mpd_vect$mpdist), 1), 24167.6)
    expect_equal(round(sum(s_mpd_vect$data[[1]]) / sd(s_mpd_vect$data[[1]]), 1), 13017.9)
    expect_equal(round(sum(s_mpd_vect$data[[2]]) / sd(s_mpd_vect$data[[2]]), 1), 1620.5)
  })

  test_that("Sub Annotation Vector", {
    s_a_vector <- a_vector[1000L:3000L]

    expect_equal(round(sum(s_a_vector$av) / sd(s_a_vector$av), 3), 8243.961)
  })

  test_that("Sub Salient", {
    s_salient <- salient[1000L:3000L]

    expect_equal(round(sum(s_salient$salient$indexes) / sd(s_salient$salient$indexes), 4), 50.5029)
    expect_equal(round(sum(s_salient$salient$idx_bit_size) / sd(s_salient$salient$idx_bit_size), 2), 95.25)
    expect_equal(sum(s_salient$salient$bits), 18)
  })

  test_that("Sub Motif", {
    s_motif <- motif[1000L:3000L]
    expect_equal(s_motif$motif$motif_idx[[1L]], c(366, 1147))
    expect_equal(s_motif$motif$motif_neighbor[[1L]], c(1348, 560, 175, 1552, 951))
  })

  test_that("Sub Discord", {
    s_discord <- discord[1000L:9000L]
    expect_equal(s_discord$discord$discord_idx[[1L]], 3891)
    expect_equal(s_discord$discord$discord_neighbor[[1L]], 633)
  })

  test_that("Head Chain", {
    h_chain <- utils::head(chain, 8000.0)
    expect_equal(sum(h_chain$chain$best), 39574.0)
    expect_equal(length(h_chain$chain$best), 6L)
  })

  test_that("Tail Chain", {
    t_chain <- utils::tail(chain, 4000)
    expect_equal(sum(t_chain$chain$best), 14625)
    expect_equal(length(t_chain$chain$best), 9)
  })

  test_that("Corrected Arc Count", {
    expect_equal(round(mean(cac$cac), 3L), 0.377)
    expect_equal(round(sd(cac$cac), 3L), 0.341)
    expect_equal(round(min(cac$cac), 3L), 0)
    expect_equal(max(cac$cac), 1.0)
  })

  test_that("Head Arc Count", {
    h_cac <- utils::head(cac, offset)
    expect_equal(round(mean(h_cac$cac), 4L), 0.7198)
    expect_equal(round(sd(h_cac$cac), 3L), 0.316)
    expect_equal(round(min(h_cac$cac), 3L), 0.209)
    expect_equal(max(h_cac$cac), 1.0)
  })

  test_that("Tail Arc Count", {
    t_cac <- utils::tail(cac, offset)
    expect_equal(round(mean(t_cac$cac), 4L), 0.7654)
    expect_equal(round(sd(t_cac$cac), 3L), 0.273)
    expect_equal(round(min(t_cac$cac), 3L), 0.255)
    expect_equal(max(t_cac$cac), 1.0)
  })

  test_that("Segments found", {
    expect_equal(segments$fluss, 4902.0)
  })

  test_that("Head Segments found", {
    h_segment <- utils::head(segments, offset)
    expect_equal(h_segment$fluss, 1228.0)
  })

  test_that("Tail Segments found", {
    t_segment <- utils::tail(segments, offset)
    expect_equal(t_segment$fluss, 1643.0)
  })
}
