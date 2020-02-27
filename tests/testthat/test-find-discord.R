if (!testthat:::on_cran()) {
  context("Testing Discord function")

  w <- 50
  data <- tsmp::mp_gait_data
  mp <- tsmp(data, window_size = w, exclusion_zone = 1 / 4, verbose = 0)
  mp <- find_discord(mp)

  test_that("Success", {
    expect_silent(find_discord(mp))
  })

  test_that("Discords", {
    expect_equal(mp$discord$discord_idx, list((48)))
    expect_equal(mp$discord$discord_neighbor, list(c(483, 184, 584)))
  })
}
