if (!identical(Sys.getenv("NOT_CRAN"), "true") ||
  identical(Sys.getenv("TRAVIS"), "true") ||
  identical(Sys.getenv("APPVEYOR"), "True")) {
  # Use minimal fonts.conf to speed up fc-cache
  gdtools::set_dummy_conf()
}
