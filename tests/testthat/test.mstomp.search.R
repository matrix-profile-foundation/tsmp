context("Testing mSTOMP Search")
library(tsmp)

w <- toy_data$sub.len
mp <- mstomp(toy_data$data[1:200, ], w, verbose = 0)
motifs <- guide.search(list(toy_data$data[1:200, 1], toy_data$data[1:200, 2], toy_data$data[1:200, 3]), w, mp$mp, mp$pi, 2)
motifs.t <- guide.search(t(toy_data$data[1:200, ]), w, mp$mp, mp$pi, 2)
motifs.u <- unconstrain.search(list(toy_data$data[1:200, 1], toy_data$data[1:200, 2], toy_data$data[1:200, 3]), w, mp$mp, mp$pi, 2)
motifs.ut <- unconstrain.search(t(toy_data$data[1:200, ]), w, mp$mp, mp$pi, 2)

test_that("Vectors", {
  expect_message(unconstrain.search(toy_data$data[1:200, 1], w, mp$mp, mp$pi, 2), regexp = "Searching")
  expect_silent(guide.search(toy_data$data[1:200, ], w, mp$mp, mp$pi, 1))
})

test_that("Errors", {
  # unknown type
  expect_error(unconstrain.search(table(toy_data$data[1:200, ]), w, mp$mp, mp$pi, 2), regexp = "Unknown type of data")
  expect_error(guide.search(table(toy_data$data[1:200, ]), w, mp$mp, mp$pi, 2), regexp = "Unknown type of data")
})

test_that("Message", {
  mpnm <- mstomp(toy_data$data[1:200, ], 100, verbose = 0)
  expect_message(unconstrain.search(toy_data$data[1:200, ], 100, mpnm$mp, mpnm$pi, 2), regexp = "No motifs found")
})

test_that("Guide Search", {
  expect_equal(motifs, motifs.t)
  expect_equal(guide.search(as.data.frame(toy_data$data[1:200, ]), w, mp$mp, mp$pi, 2), motifs)
  expect_equal(motifs$motif.idx, c(44, 108))
  expect_equal(motifs$motif.dim, list(c(2, 3), c(2, 3)))
})

test_that("Unguide Search", {
  expect_equal(motifs.u, motifs.ut)
  expect_equal(unconstrain.search(as.data.frame(toy_data$data[1:200, ]), w, mp$mp, mp$pi, 2), motifs.u)
  expect_equal(motifs.u$motif.idx, c(33, 102, 57, 9, 127, 81))
  expect_equal(motifs.u$motif.dim, list(1, 1, 1, 1, 1, 3))
})
