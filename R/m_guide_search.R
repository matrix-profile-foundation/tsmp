#' Guided MOTIF Discovery for Multidimensional Matrix Profile
#'
#' @details
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
#' w <- mp_toy_data$sub.len
#' mp <- mstomp(mp_toy_data$data[1:200,], w, verbose = 0)
#' motifs <- guide.search(mp_toy_data$data[1:200,], w, mp$mp, mp$pi, 2)
#' \dontrun{
#' w <- mp_toy_data$sub.len
#' mp <- mstomp.par(mp_toy_data$data, w, verbose = 0)
#' motifs <- guide.search(mp_toy_data$data, w, mp$mp, mp$pi, 2)
#' }

guide_search <- function(data, window_size, matrix_profile, profile_index, n_dim) {

  ## transform data list into matrix
  if (is.matrix(data) || is.data.frame(data)) {
    if (is.data.frame(data)) {
      data <- as.matrix(data)
    } # just to be uniform
    if (ncol(data) > nrow(data)) {
      data <- t(data)
    }
  } else if (is.list(data)) {
    data_len <- length(data[[1]])
    data_dim <- length(data)

    for (i in 1:data_dim) {
      len <- length(data[[i]])
      # Fix TS size with NaN
      if (len < data_len) {
        data[[i]] <- c(data[[i]], rep(NA, data_len - len))
      }
    }
    # transform data into matrix (each column is a TS)
    data <- sapply(data, cbind)
  } else if (is.vector(data)) {
    # transform data into 1-col matrix
    data <- as.matrix(data) # just to be uniform
  } else {
    stop("Error: Unknown type of data. Must be: matrix, data_frame, vector or list.", call. = FALSE)
  }

  matrix_profile <- matrix_profile[, n_dim]
  profile_index <- profile_index[, n_dim]
  motif_idx <- which.min(matrix_profile)
  motif_idx <- sort(c(motif_idx, profile_index[motif_idx]))

  motif_1 <- as.matrix(data[motif_idx[1]:(motif_idx[1] + window_size - 1), ]) # as.matrix(): hack for vectors
  motif_2 <- as.matrix(data[motif_idx[2]:(motif_idx[2] + window_size - 1), ]) # as.matrix(): hack for vectors

  motif_dim <- sort(apply(abs(motif_1 - motif_2), 2, sum), index.return = TRUE)$ix
  motif_dim <- sort(motif_dim[1:n_dim])
  motif_dim <- list(motif_dim, motif_dim)

  return(list(motif_idx = motif_idx, motif_dim = motif_dim))
}
