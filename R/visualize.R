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
#' result <- compute(mp_toy_data$data[, 1], 80)
#' visualize(result)
visualize <- function(obj) {
  if (!is.null(obj$profile)) {
    invisible(graphics::plot(obj$profile))
  } else {
    invisible(graphics::plot(obj))
  }
}
