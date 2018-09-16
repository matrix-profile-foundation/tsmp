if (identical(Sys.getenv("NOT_CRAN"), "true") &&
  !identical(Sys.getenv("TRAVIS"), "true") &&
  !identical(Sys.getenv("APPVEYOR"), "True")) {
  if (requireNamespace("spelling", quietly = TRUE)) {
    spelling::spell_check_test(vignettes = TRUE, error = FALSE)
  }
}
