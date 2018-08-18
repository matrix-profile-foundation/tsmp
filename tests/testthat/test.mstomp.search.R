context("Testing FLUSS functions")
library(tsmp)

w <- toy_data$sub.len
mp <- mstomp(toy_data$data[1:200, ], w, verbose = 0)
motifs <- guide.search(toy_data$data[1:200, ], w, mp$mp, mp$pi, 2)
motifs.u <- unconstrain.search(toy_data$data[1:200, ], w, mp$mp, mp$pi, 2)


test_that("Guide Search", {
  expect_equal(motifs$motif.idx, c(44, 108))
  expect_equal(motifs$motif.dim, list(c(2,3), c(2,3)))
})

test_that("Unguide Search", {
  expect_equal(motifs.u$motif.idx, c(33, 102, 57, 9, 127, 81))
  expect_equal(motifs.u$motif.dim, list(1, 1, 1, 1, 1, 3))
})
