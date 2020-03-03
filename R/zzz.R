tsmp_default_options <- list(
  tsmp.verbose = 2,
  tsmp.exclusion_zone = 1 / 2
)

.onLoad <- function(libname, pkgname) {
  op <- options()
  toset <- !(names(tsmp_default_options) %in% names(op))
  if (any(toset)) {
    options(tsmp_default_options[toset])
  }

  invisible()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to Matrix Profile")
}
