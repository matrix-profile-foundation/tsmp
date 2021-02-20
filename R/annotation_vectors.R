#' Computes the annotation vector that favors number of zero crossing
#'
#' @param .mp a Matrix Profile object.
#' @param data a `vector` or a column `matrix` of `numeric`.
#' @param apply logical. (Default is `FALSE`). Applies the Annotation Vector over the Matrix Profile.
#'  Use with caution.
#'
#' @return Returns the input `.mp` object with an embedded annotation vector.
#' @export
#' @family Annotation vectors
#' @references * Dau HA, Keogh E. Matrix Profile V: A Generic Technique to Incorporate Domain
#'   Knowledge into Motif Discovery. In: Proceedings of the 23rd ACM SIGKDD International Conference
#'   on Knowledge Discovery and Data Mining - KDD '17. New York, New York, USA: ACM Press; 2017. p.
#'   125-34.
#' @examples
#' data <- mp_test_data$train$data[1:1000]
#' w <- 50
#' mp <- tsmp(data, window_size = w, verbose = 0)
#' av <- av_zerocrossing(mp, apply = TRUE)
av_zerocrossing <- function(.mp, data, apply = FALSE) {
  if ("Valmod" %in% class(.mp)) {
    stop("Function not implemented for objects of class `Valmod`.")
  }

  # validate input ----
  checkmate::qassert(.mp, "L+")
  if (missing(data) && !is.null(.mp$data)) {
    data <- as.numeric(.mp$data[[1L]])
    checkmate::qassert(data, "N+")
  }
  window_size <- .mp$w

  av <- matrixprofiler::zero_crossing(data, window_size)
  .mp$av <- matrixprofiler::normalize(av, 0, 1)

  class(.mp) <- update_class(class(.mp), "AnnotationVector")

  if (apply == TRUE) {
    .mp <- av_apply(.mp)
  }

  return(.mp)
}
