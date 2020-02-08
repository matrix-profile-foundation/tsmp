#' Plots an object generated from one of the algorithms. In some cases multiple plots will be generated
#'
#' @param profile a `MatrixProfile` or `PMP` object.
#'
#' @export
#'
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @family Main API
#'
#' @examples
#' result <- compute(mp_toy_data$data[, 1], 80)
#' visualize(result)
visualize <- function(profile) {
  if (!is.null(profile$profile)) {
    invisible(graphics::plot(profile$profile))
  } else {
    invisible(graphics::plot(profile))
  }
}
