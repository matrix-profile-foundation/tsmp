#' Moving window minimum
#'
#' @param data
#' @param window_size
#'
#' @return
#' @keywords internal
#' @noRd
movmin <- function(data, window_size) {
  data <- as.vector(data)
  data_size <- length(data)
  window_size <- as.integer(window_size)
  if (window_size <= 1) {
    return(data)
  }
  if (window_size > data_size) window_size <- data_size

  y <- .C("movmin", as.double(data),
    y = double(data_size), as.integer(data_size), as.integer(window_size),
    NAOK = TRUE, PACKAGE = "tsmp"
  )$y

  return(y)
}

#' Moving window maximum
#'
#' @param data
#' @param window_size
#'
#' @return
#' @keywords internal
#' @noRd
movmax <- function(data, window_size) {
  data <- as.vector(data)
  data_size <- length(data)
  window_size <- as.integer(window_size)
  if (window_size <= 1) {
    return(data)
  }
  if (window_size > data_size) window_size <- data_size
  y <- double(data_size)

  y <- .C("movmax", as.double(data),
    y = double(data_size), as.integer(data_size), as.integer(window_size),
    NAOK = TRUE, PACKAGE = "tsmp"
  )$y

  return(y)
}
