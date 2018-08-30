#' Guided MOTIF Discovery for Multidimensional Matrix Profile
#'
#' Guided MOTIF Discovery for Multidimensional Matrix Profile
#'
#' Although this functions handles Multivariate Time Series, it can also be used to handle
#' Univariate Time Series.
#'
#' @param data a `matrix` of `numeric`, where each column is a time series. Accepts `vector` (see
#'   details), `list` and `data.frame` too.
#' @param window.size an `int` with the size of the sliding window.
#' @param matrix.profile multidimensional matrix profile (from [mstomp()] or [mstomp.par()]).
#' @param profile.index multidimensional profile index (from [mstomp()] or [mstomp.par()]).
#' @param n.dim an `int`. The dimensionality of the MOTIF to find.
#'
#' @return Returns the `motif.idx` with the index of MOTIFs founded and `motif.dim` with the spanned
#'   dimensions of respective MOTIF.
#'
#' @export
#'
#' @family mstomp
#' @seealso [mstomp()], [mstomp.par()], [unconstrain.search()]
#' @references * Yeh CM, Kavantzas N, Keogh E. Matrix Profile VI : Meaningful Multidimensional Motif
#'   Discovery.
#' @references * Zhu Y, Imamura M, Nikovski D, Keogh E. Matrix Profile VII: Time Series Chains: A
#'   New Primitive for Time Series Data Mining. Knowl Inf Syst. 2018 Jun 2;1–27.
#' @references Website: <https://sites.google.com/view/mstamp/>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#' # This is a fast toy example and results are useless. For a complete result, run the code inside
#' #'Not run' section below.
#' w <- toy_data$sub.len
#' mp <- mstomp(toy_data$data[1:200,], w, verbose = 0)
#' motifs <- guide.search(toy_data$data[1:200,], w, mp$mp, mp$pi, 2)
#' \dontrun{
#' w <- toy_data$sub.len
#' mp <- mstomp.par(toy_data$data, w, verbose = 0)
#' motifs <- guide.search(toy_data$data, w, mp$mp, mp$pi, 2)
#' }

guide.search <- function(data, window.size, matrix.profile, profile.index, n.dim) {

  ## transform data list into matrix
  if (is.matrix(data) || is.data.frame(data)) {
    if (is.data.frame(data)) {
      data <- as.matrix(data)
    } # just to be uniform
    if (ncol(data) > nrow(data)) {
      data <- t(data)
    }
  } else if (is.list(data)) {
    data.len <- length(data[[1]])
    data.dim <- length(data)

    for (i in 1:data.dim) {
      len <- length(data[[i]])
      # Fix TS size with NaN
      if (len < data.len) {
        data[[i]] <- c(data[[i]], rep(NA, data.len - len))
      }
    }
    # transform data into matrix (each column is a TS)
    data <- sapply(data, cbind)
  } else if (is.vector(data)) {
    # transform data into 1-col matrix
    data <- as.matrix(data) # just to be uniform
  } else {
    stop("Error: Unknown type of data. Must be: matrix, data.frame, vector or list")
  }

  matrix.profile <- matrix.profile[, n.dim]
  profile.index <- profile.index[, n.dim]
  motif.idx <- which.min(matrix.profile)
  motif.idx <- sort(c(motif.idx, profile.index[motif.idx]))

  motif.1 <- as.matrix(data[motif.idx[1]:(motif.idx[1] + window.size - 1), ]) # as.matrix(): hack for vectors
  motif.2 <- as.matrix(data[motif.idx[2]:(motif.idx[2] + window.size - 1), ]) # as.matrix(): hack for vectors

  motif.dim <- sort(apply(abs(motif.1 - motif.2), 2, sum), index.return = TRUE)$ix
  motif.dim <- sort(motif.dim[1:n.dim])
  motif.dim <- list(motif.dim, motif.dim)

  return(list(motif.idx = motif.idx, motif.dim = motif.dim))
}
