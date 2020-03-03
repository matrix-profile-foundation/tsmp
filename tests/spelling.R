if (requireNamespace("spelling", quietly = TRUE)) {
  if (!covr::in_covr()) {
    spelling::spell_check_test(
      vignettes = TRUE, error = FALSE,
      skip_on_cran = TRUE
    )
  }
}
