library(testthat)

if (all(skip_on_cran(), skip_on_travis(), (is.null(skip_on_appveyor())))) {
  if (requireNamespace("spelling", quietly = TRUE)) {
    spelling::spell_check_test(vignettes = TRUE, error = FALSE)
  }
}
