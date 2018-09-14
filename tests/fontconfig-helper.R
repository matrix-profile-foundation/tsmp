# Use minimal fonts.conf to speed up fc-cache
if (!skip_on_cran() || !skip_on_travis() || !(is.null(skip_on_appveyor()))) {
  gdtools::set_dummy_conf()
}
