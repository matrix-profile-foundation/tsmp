skip_on_cran()
skip_if_not_installed("vdiffr")
library(tsmp)

data <- mp_fluss_data$tilt_abp$data[1:1000]
w <- 10
truth <- 400
nseg <- 3
mp <- tsmp(data, window_size = w, verbose = 0)
cac <- fluss_cac(mp)
segments <- fluss_extract(cac, nseg)
chain <- find_chains(mp)
motif <- find_motif(mp)
discord <- find_discord(mp)
mps <- salient_subsequences(mp, n_bits = c(4, 6, 8), verbose = 0)
val <- tsmp(data, window_size = c(30, 40), verbose = 0, mode = "valmod")
val_motif <- find_motif(val)
snippet <- find_snippet(mp_gait_data[1:500, ], 80)

plot_arcs_test <- function() plot_arcs(pairs = matrix(c(5, 10, 1, 10, 20, 5), ncol = 2, byrow = TRUE))
plot_arccount_test <- function() plot.ArcCount(cac)
plot_matrixprofile_test <- function() plot.MatrixProfile(mp)
plot_fluss_test <- function() plot.Fluss(segments)
plot_chain_test <- function() plot.Chain(chain)
plot_motif_test <- function() plot.Motif(motif)
plot_discord_test <- function() plot.Discord(discord)
# cmdscale() may reflect either axis across platforms. Orient each axis using
# its largest absolute coordinate; the signs retain the existing SVG reference.
# This affects only this visual test, not salient_mds() or plot.Salient().
orient_salient_mds <- function(coords) {
  target_sign <- c(1, -1)
  for (axis in seq_len(ncol(coords))) {
    anchor <- which.max(abs(coords[, axis]))
    if (coords[anchor, axis] * target_sign[axis] < 0) {
      coords[, axis] <- -coords[, axis]
    }
  }
  coords
}

plot_salient_test <- function() {
  plot <- plot.Salient
  environment(plot) <- list2env(list(
    salient_mds = function(...) orient_salient_mds(salient_mds(...))
  ), parent = environment(plot))
  plot(mps)
}
plot_valmod_test <- function() plot.Valmod(val)
plot_snippet_test <- function() plot.Snippet(snippet)

mdata <- mp_toy_data$data[1:200, ]
mw <- mp_toy_data$sub_len
mmp <- tsmp(mdata, window_size = mw, mode = "mstomp", verbose = 0)
smp <- tsmp(mdata, window_size = mw, mode = "simple", verbose = 0)
mmotif <- find_motif(mmp, n_motifs = 2)

plot_multimatrixprofile_test <- function() plot.MultiMatrixProfile(mmp)
plot_multimotif_test <- function() plot.MultiMotif(mmotif)
plot_simplematrixprofile_test <- function() plot.SimpleMatrixProfile(smp)

test_that("Plot", {
  vdiffr::expect_doppelganger("plot arcs", plot_arcs_test)
  vdiffr::expect_doppelganger("plot arc count", plot_arccount_test)
  vdiffr::expect_doppelganger("plot matrix profile", plot_matrixprofile_test)
  vdiffr::expect_doppelganger("plot fluss", plot_fluss_test)
  vdiffr::expect_doppelganger("plot chain", plot_chain_test)
  vdiffr::expect_doppelganger("plot motif", plot_motif_test)
  vdiffr::expect_doppelganger("plot discord", plot_discord_test)
  vdiffr::expect_doppelganger("plot multi matrix profile", plot_multimatrixprofile_test)
  vdiffr::expect_doppelganger("plot multimotif", plot_multimotif_test)
  vdiffr::expect_doppelganger("plot simple matrix profile", plot_simplematrixprofile_test)
  vdiffr::expect_doppelganger("plot salient", plot_salient_test)
  vdiffr::expect_doppelganger("plot valmod", plot_valmod_test)
  vdiffr::expect_doppelganger("plot Snippet", plot_snippet_test)
})

test_that("Print", {
  expect_snapshot_output(print(segments))
  expect_snapshot_output(print(cac))
  expect_snapshot_output(print(mp))
  expect_snapshot_output(print(mmp))
  expect_snapshot_output(print(smp))
  expect_snapshot_output(print(chain))
  expect_snapshot_output(print(motif))
  expect_snapshot_output(print(discord))
  expect_snapshot_output(print(mmotif))
  expect_snapshot_output(print(mps))
  expect_snapshot_output(print(val))
  expect_snapshot_output(print(val_motif))
  expect_snapshot_output(print(snippet))
})


test_that("Salient plot orientation ignores reflections but preserves geometry", {
  coords <- salient_mds(mps)
  oriented <- orient_salient_mds(coords)
  for (x_sign in c(-1, 1)) {
    for (y_sign in c(-1, 1)) {
      reflected <- sweep(coords, 2L, c(x_sign, y_sign), "*")
      expect_equal(orient_salient_mds(reflected), oriented)
    }
  }
  expect_equal(as.matrix(dist(oriented)), as.matrix(dist(coords)))
  # A change in scale must remain visible to the snapshot comparison.
  expect_equal(orient_salient_mds(2 * coords), 2 * oriented)
})
