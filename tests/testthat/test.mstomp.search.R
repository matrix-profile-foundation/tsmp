context("Testing mSTOMP Search")
library(tsmp)

w <- toy_data$sub_len
mp <- mstomp(toy_data$data[1:200, ], w, verbose = 0)
motifs <- guide_search(list(toy_data$data[1:200, 1], toy_data$data[1:200, 2], toy_data$data[1:200, 3]), w, mp$mp, mp$pi, 2)
motifs_t <- guide_search(t(toy_data$data[1:200, ]), w, mp$mp, mp$pi, 2)
motifs_u <- unconstrain_search(list(toy_data$data[1:200, 1], toy_data$data[1:200, 2], toy_data$data[1:200, 3]), w, mp$mp, mp$pi, 2)
motifs_ut <- unconstrain_search(t(toy_data$data[1:200, ]), w, mp$mp, mp$pi, 2)

test_that("Vectors", {
  expect_message(unconstrain_search(toy_data$data[1:200, 1], w, mp$mp, mp$pi, 2), regexp = "Searching")
  expect_silent(guide_search(toy_data$data[1:200, ], w, mp$mp, mp$pi, 1))
})

test_that("Errors", {
  # unknown type
  expect_error(unconstrain_search(table(toy_data$data[1:200, ]), w, mp$mp, mp$pi, 2), regexp = "Unknown type of data")
  expect_error(guide_search(table(toy_data$data[1:200, ]), w, mp$mp, mp$pi, 2), regexp = "Unknown type of data")
})

test_that("Message", {
  mpnm <- mstomp(toy_data$data[1:200, ], 100, verbose = 0)
  expect_message(unconstrain_search(toy_data$data[1:200, ], 100, mpnm$mp, mpnm$pi, 2), regexp = "No motifs found")
})

test_that("Guide Search", {
  expect_equal(motifs, motifs_t)
  expect_equal(guide_search(as.data.frame(toy_data$data[1:200, ]), w, mp$mp, mp$pi, 2), motifs)
  expect_equal(motifs$motif_idx, c(44, 108))
  expect_equal(motifs$motif_dim, list(c(2, 3), c(2, 3)))
})

test_that("Unguide Search", {
  expect_equal(motifs_u, motifs_ut)
  expect_equal(unconstrain_search(as.data.frame(toy_data$data[1:200, ]), w, mp$mp, mp$pi, 2), motifs_u)
  expect_equal(motifs_u$motif_idx, c(33, 102, 57, 9, 127, 81))
  expect_equal(motifs_u$motif_dim, list(1, 1, 1, 1, 1, 3))
})
