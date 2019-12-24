if (!testthat:::on_cran()) {
  library(vdiffr)
  library(tsmp)

  context("Testing Plots")

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
  plot_salient_test <- function() plot.Salient(mps)
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
    expect_doppelganger("plot arcs", plot_arcs_test)
    expect_doppelganger("plot arc count", plot_arccount_test)
    expect_doppelganger("plot matrix profile", plot_matrixprofile_test)
    expect_doppelganger("plot fluss", plot_fluss_test)
    expect_doppelganger("plot chain", plot_chain_test)
    expect_doppelganger("plot motif", plot_motif_test)
    expect_doppelganger("plot discord", plot_discord_test)
    expect_doppelganger("plot multi matrix profile", plot_multimatrixprofile_test)
    expect_doppelganger("plot multimotif", plot_multimotif_test)
    expect_doppelganger("plot simple matrix profile", plot_simplematrixprofile_test)
    expect_doppelganger("plot salient", plot_salient_test)
    expect_doppelganger("plot valmod", plot_valmod_test)
    expect_doppelganger("plot Snippet", plot_snippet_test)
  })

  context("Testing Print")
  upd <- FALSE
  if (is_testing()) {
    path <- "../prints/"
  } else {
    path <- "./tests/prints/"
  }

  test_that("Print", {
    expect_known_output(segments, file = paste0(path, "fluss-print"), print = TRUE, update = upd)
    expect_known_output(cac, file = paste0(path, "cac-print"), print = TRUE, update = upd)
    expect_known_output(mp, file = paste0(path, "mp-print"), print = TRUE, update = upd)
    expect_known_output(mmp, file = paste0(path, "mmp-print"), print = TRUE, update = upd)
    expect_known_output(smp, file = paste0(path, "smp-print"), print = TRUE, update = upd)
    expect_known_output(chain, file = paste0(path, "chain-print"), print = TRUE, update = upd)
    expect_known_output(motif, file = paste0(path, "motif-print"), print = TRUE, update = upd)
    expect_known_output(discord, file = paste0(path, "discord-print"), print = TRUE, update = upd)
    expect_known_output(mmotif, file = paste0(path, "mmotif-print"), print = TRUE, update = upd)
    expect_known_output(mps, file = paste0(path, "salient-print"), print = TRUE, update = upd)
    expect_known_output(val, file = paste0(path, "valmod-print"), print = TRUE, update = upd)
    expect_known_output(val_motif, file = paste0(path, "motif_valmod-print"), print = TRUE, update = upd)
    expect_known_output(snippet, file = paste0(path, "snippet-print"), print = TRUE, update = upd)
  })
}
