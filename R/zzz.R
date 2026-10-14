.onLoad <- function(libname, pkgname) { # nolint
  tryCatch(debugme::debugme(), error = identity)

  invisible()
}

.onAttach <- function(libname, pkgname) { # nolint
  packageStartupMessage("Welcome to Matrix Profile")
}

.onUnload <- function(libpath) { # nolint
  unloadNamespace("debugme")
}

#' Returns a message to the user about function deprecation.
#'
#' Adapted from gg_dep from package ggplot2
#'
#' @param version The last version of matrixprofiler where this function was good (in other words,
#' the last version where it was not deprecated).
#' @param msg The message to print.
#'
#' @return Invisibly returns nothing
#' @keywords internal
#' @noRd
mp_dep <- function(version, msg) {
  v <- as.package_version(version)
  cv <- utils::packageVersion("matrixprofiler")

  # If current major number is greater than last-good major number, or if
  # current minor number is more than 1 greater than last-good minor number,
  # return an error.
  if (cv[[1, 1]] > v[[1, 1]] || cv[[1, 2]] > v[[1, 2]] + 1) {
    stop(msg, " (Defunct; last used in version ", version, ")",
      call. = FALSE
    )

    # If minor number differs by one, give a warning
  } else if (cv[[1, 2]] > v[[1, 2]]) {
    warning(msg, " (Deprecated; last used in version ", version, ")",
      call. = FALSE
    )

    # If only subminor number is greater, provide a message
  } else if (cv[[1, 3]] > v[[1, 3]]) {
    message(msg, " (Deprecated; last used in version ", version, ")")
  }

  invisible()
}
