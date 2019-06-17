#' Time Series Snippets: A New Primitive for Time Series Data Mining
#'
#' MPdist is a recently introduced distance measure which considers two time series to be similar
#'
#' @details MPdist
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
#' snippet within time sereis
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
find_snippet <- function(data, s_size, n_snippets = 2, window_size = s_size / 2) {

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

  if (s_size < 4) {
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

  ## initialization
  distances <- NULL
  indexes <- NULL
  snippet <- NULL
  minis <- Inf
  fraction <- NULL
  snippetidx <- NULL
  distances_snippet <- NULL

  len <- s_size * ceiling(nrow(data) / s_size) - nrow(data)
  ## padding zeros to the end of the time series.
  data <- rbind(data, matrix(0, nrow = len))

  ## compute all the profiles
  for (i in seq.int(1, nrow(data) - s_size, s_size)) {
    indexes <- c(indexes, i)
    distance <- mpdist(data, data[i:(i + s_size - 1), , drop = FALSE], window_size, type = "vector")
    distances <- rbind(distances, distance$mpdist)
  }

  ## calculating the Nth snippets
  for (n in 1:n_snippets) {
    minims <- Inf
    for (i in seq_len(nrow(distances))) {
      if (minims > sum(pmin(distances[i, ], minis))) {
        minims <- sum(pmin(distances[i, ], minis))
        index <- i
      }
    }
    minis <- pmin(distances[index, ], minis)
    snippet <- c(snippet, data[indexes[index]:(indexes[index] + s_size - 1)])
    snippetidx <- c(snippetidx, indexes[index])

    distance <- mpdist(data, data[indexes[index]:(indexes[index] + s_size - 1), , drop = FALSE], window_size, type = "vector")
    distances_snippet <- rbind(distances_snippet, distance$mpdist)

    graphics::plot(distance$mpdist,
      type = "l", ylab = "distance",
      main = paste0("MPdist-", n, " location-", indexes[index])
    )

    graphics::plot(data[indexes[index]:(indexes[index] + s_size - 1)],
      type = "l", ylab = "data",
      main = paste0("snippet-", n, " location-", indexes[index])
    )
  }
  ## Calculating the fraction of each snippet
  totalmin <- do.call(pmin, as.data.frame(t(distances_snippet)))

  for (i in 1:n_snippets) {
    a <- (distances_snippet[i, ] <= totalmin)
    fraction <- c(fraction, sum(a) / (nrow(data) - s_size))
    # # # #     just in case we have the same value for both
    totalmin[a] <- totalmin[a] - 1
  }

  ## Calculating the horizontal regime bar for 2 snippets
  if (n_snippets == 2) {
    totalmin <- do.call(pmin, as.data.frame(t(distances_snippet)))

    a <- (distances_snippet[1, ] <= totalmin)

    for (i in seq.int(1, length(a) - s_size, s_size)) {
      if (sum(a[i:(i + s_size - 1)]) > (1 / 2 * s_size)) {
        a[i:(i + s_size - 1)] <- TRUE
      } else {
        a[i:(i + s_size - 1)] <- FALSE
      }
    }

    # # # # # # # # # # # # # # # # #
    c <- matrix(0, nrow = length(a), ncol = 3)
    j <- 1
    i <- 1

    while (i <= length(a)) {
      b <- which(a[(i + 1):length(a)] != a[i])[1]
      if (is.na(b)) {
        break
      } else {
        c[j, ] <- c(i, i + b - 1, a[i])
        j <- j + 1
        i <- i + b
      }
    }
    if (i <= length(a)) {
      c[j, ] <- c(i, length(a), a[i])
    }
    c <- c[1:j, , drop = FALSE]

    graphics::plot(0.5, 0.5,
      ylab = "", xlab = "Index",
      type = "n", main = "Horizontal regime bar", xlim = c(min(c[, 1]), max(c[, 2])), ylim = c(0, 2)
    )

    for (i in 1:nrow(c)) {
      d <- (c[i, 1]:c[i, 2])
      if (c[i, 3] == 0) {
        graphics::lines(d, rep(1, length(d)), col = "blue")
      } else {
        graphics::lines(d, rep(1, length(d)), col = "red")
      }
    }
  }

  return(list(snippets = snippet, snippetsidx = snippetidx, fractions = fraction))
}
