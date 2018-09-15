library(testthat)
if (!skip_on_cran() || !skip_on_travis() || !(is.null(skip_on_appveyor()))) {
  # Use minimal fonts.conf to speed up fc-cache
  gdtools::set_dummy_conf()
}
