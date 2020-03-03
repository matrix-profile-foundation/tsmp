#' Time Series Snippets: A New Primitive for Time Series Data Mining
#'
#' Time Series Snippets tries to solve mainly the common problem of summarization "Show me some
#' representative/typical data". As stated by the original paper, potential uses of snippets are:
#' integrating summarizations of files directly into an operating, production of automatically
#' generated reports, for example, summarize a sleep study and also can be used to support a
#' host of higher-level tasks, including the comparison of massive data collections.
#'
#' Motifs vs. snippets: While motifs reward fidelity of conservation, snippets also rewards coverage.
#' Informally, coverage is some measure of how much of the data is explained or represented by
#' a given snippet.
#'
#' Shapelets vs. snippets: shapelets are defined as subsequences that are maximally representative
#' of a class. Shapelets are supervised, snippets are unsupervised. Shapelets are generally biased
#' to be as short as possible. In contrast, we want snippets to be longer, to intuitively capture
#' the "flavor" of the time series.
#'
#' @param data a `matrix` or a `vector`.
#' @param s_size an int. Size of snippet.
#' @param n_snippets an `int`. Number of snippets to find. (Default is `2`).
#' @param window_size an `int`. The size of the sliding window used to compare the data. Must be smaller
#' than `s_size`. (Default is `s_size / 2`).
#'
#' @return Returns the snippet : a list of n_snippets snippets
#' fraction : fraction of each snippet
#' snippetidx : the location of each
#' snippet within time series
#'
#' @export
#'
#' @references * Imani S, Madrid F, Ding W, Crouter S, Keogh E. Matrix Profile XIII:
#' Time Series Snippets: A New Primitive for Time Series Data Mining. In: 2018 IEEE International
#' Conference on Data Mining (ICDM). 2018.
#' @references * Gharghabi S, Imani S, Bagnall A, Darvishzadeh A, Keogh E. Matrix Profile XII:
#' MPdist: A Novel Time Series Distance Measure to Allow Data Mining in More Challenging Scenarios.
#' In: 2018 IEEE International Conference on Data Mining (ICDM). 2018.
#' @references Website: <https://sites.google.com/site/snippetfinder/>
#'
#' @examples
#'
#' snippets <- find_snippet(mp_fluss_data$walkjogrun$data[1:300], 40, n_snippets = 2)
#' \dontrun{
#' snippets <- find_snippet(mp_fluss_data$walkjogrun$data, 120, n_snippets = 3)
#' plot(snippets)
#' }
#'
find_snippet <- function(data, s_size, n_snippets = 2L, window_size = s_size / 2L) {

  # currently is about 3x slower than MATLAB. Not bad for R.
  # transform data into matrix
  if (is.vector(data)) {
    data <- as.matrix(data)
  }
  else if (is.matrix(data)) {
    if (ncol(data) > nrow(data)) {
      data <- t(data)
    }
  } else {
    stop("Unknown type of data. Must be: a column matrix or a vector.", call. = FALSE)
  }

  if (s_size < 4.0) {
    stop("`s_size` must be at least 4.", call. = FALSE)
  }

  ## check input
  if (nrow(data) < (2 * s_size)) {
    stop("Error: Time series is too short relative to desired snippet Length", call. = FALSE)
  }

  window_size <- floor(window_size)

  if (window_size >= s_size) {
    stop("Error: `window_size` must be smaller than `s_size`.", call. = FALSE)
  }

  #### padding zeros to the end of the time series. ####
  pad <- s_size * ceiling(nrow(data) / s_size) - nrow(data)
  data <- rbind(data, matrix(0, nrow = pad))

  #### compute all the profiles ####
  indexes <- seq.int(1, nrow(data) - s_size, s_size)
  distances <- matrix(nrow = length(indexes), ncol = (nrow(data) - s_size + 1))
  distances_snippet <- matrix(nrow = n_snippets, ncol = (nrow(data) - s_size + 1))

  for (j in seq_along(indexes)) {
    i <- indexes[j]
    distance <- mpdist(data, data[i:(i + s_size - 1), , drop = FALSE], window_size, type = "vector")
    distances[j, ] <- distance$mpdist
  }

  #### calculating the Nth snippets ####
  minis <- Inf
  snippetidx <- NULL
  for (n in 1:n_snippets) {
    minims <- Inf
    index <- NULL
    for (i in seq_len(nrow(distances))) {
      s <- sum(pmin(distances[i, ], minis)) # area under the profile (maximize coverage)
      if (minims > s) {
        minims <- s
        index <- i
      }
    }
    minis <- pmin(distances[index, ], minis)
    snippetidx <- c(snippetidx, indexes[index])
    distances_snippet[n, ] <- distances[index, ]
  }

  #### Calculating the fraction of each snippet ####
  totalmin <- do.call(pmin, as.data.frame(t(distances_snippet))) # col wise
  horizontal <- rep(0, length(totalmin))
  fraction <- NULL

  for (i in 1:n_snippets) {
    a <- (distances_snippet[i, ] <= totalmin)
    fraction <- c(fraction, sum(a) / (nrow(data) - s_size + 1))

    # just in case we have the same value for both
    totalmin[a] <- totalmin[a] - 1

    for (j in indexes) {
      if (sum(a[j:(j + s_size - 1)]) > (1 / 2 * s_size)) {
        a[j:(j + s_size - 1)] <- TRUE
      } else {
        a[j:(j + s_size - 1)] <- FALSE
      }
    }

    horizontal[a] <- i
  }

  ## assert percentages
  if (round(sum(fraction), 3) != 1) {
    message("DEBUG: ", round(sum(fraction), 3))
  }

  obj <- list(snippet_idx = snippetidx, snippet_frac = fraction, snippet_size = s_size, regime = horizontal, data = list(data))
  class(obj) <- "Snippet"

  return(obj)
}
