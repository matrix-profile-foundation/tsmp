if (!testthat:::on_cran()) {
  context("Testing SiMPle Fast")
  library(tsmp)

  w <- 30
  data <- mp_toy_data$data[1:250, ] # 3 dimensions matrix
  query <- mp_toy_data$data[251:500, ] # 3 dimensions matrix

  test_that("Errors", {
    # big window size
    expect_error(tsmp(data, window_size = 500, mode = "simple"), "too short relative")

    # short window size
    expect_error(simple_fast(data, window_size = 2), "must be at least 4")

    # invalid window
    expect_error(simple_fast(data, window_size = data), "window_size")

    # data and query dim must be the same
    expect_error(simple_fast(data, data[, 1], window_size = w), "Data and query dimensions")
    expect_error(simple_fast(data[, 1], data, window_size = w), "Data and query dimensions")

    # Unknown type
    expect_error(simple_fast(table(data), window_size = w), "Unknown type of data")
    expect_error(simple_fast(data, table(data), window_size = w), "Unknown type of query")
  })

  test_that("Messages", {
    expect_message(simple_fast(data, window_size = w, verbose = 1), "Finished")
    expect_warning(tsmp(data, data, data, window_size = w, verbose = 0, mode = "simple"), "Only the first two")
  })

  result_self <- simple_fast(data, window_size = w, verbose = 2)
  result_join <- simple_fast(data, query, window_size = w, verbose = 2)

  test_that("SiMPle Results", {
    expect_equal(result_self$n_dim, 3)
    expect_equal(round(sum(result_self$mp), 3), 419.509)
    expect_equal(round(sd(result_self$mp), 3), 0.841)
    expect_equal(sum(result_self$pi), 23878)
    expect_equal(round(sd(result_self$pi), 3), 64.977)

    expect_equal(result_join$n_dim, 3)
    expect_equal(round(sum(result_join$mp), 3), 908.248)
    expect_equal(round(sd(result_join$mp), 3), 2.263)
    expect_equal(sum(result_join$pi), 24981)
    expect_equal(round(sd(result_join$pi), 3), 61.021)
  })

  test_that("Find Motif", {
    result <- find_motif(result_self, data)
    expect_equal(result$motif$motif_idx, list(c(73, 89), c(160, 181), c(1, 16)))
    expect_equal(result$motif$motif_neighbor, list(numeric(), numeric(), c(120, 183, 594, 389)))
  })

  test_that("Find Discord", {
    result <- find_discord(result_self, data)
    expect_equal(result$discord$discord_idx, list((142)))
    expect_equal(result$discord$discord_neighbor, list(c(414, 378, 85)))
  })
}
