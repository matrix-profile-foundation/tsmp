#' Plots an object generated from one of the algorithms. In some cases multiple plots will be generated
#'
#' @param obj  dict_like The object to plot.
#'
#' @return
#' A list of matplotlib figures
#'
#' @export
#'
#' @family Main API
#'
#' @examples
#'
visualize <- function(obj) {
  if (!is.null(obj$profile)) {
    invisible(plot(obj$profile))
  } else {
    invisible(plot(obj))
  }
}
