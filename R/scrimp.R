#' Anytime univariate SCRIMP algorithm (experimental)
#'
#' Computes the best so far Matrix Profile and Profile Index for Univariate Time Series.
#' DISCLAIMER: This algorithm still in development by its authors.
#' Join similarity, RMP and LMP not implemented yet.
#'
#' @details
#' The Matrix Profile, has the potential to revolutionize time series data mining because of its
#' generality, versatility, simplicity and scalability. In particular it has implications for time
#' series motif discovery, time series joins, shapelet discovery (classification), density
#' estimation, semantic segmentation, visualization, rule discovery, clustering etc. The anytime
#' SCRIMP computes the Matrix Profile and Profile Index in such manner that it can be stopped before
#' its complete calculation and return the best so far results allowing ultra-fast approximate
#' solutions. `verbose` changes how much information is printed by this function; `0` means nothing,
#' `1` means text, `2` means text and sound. `exclusion.zone` is used to avoid  trivial matches.
#'
#' @param ... a `matrix` or a `vector`.
#' @param window.size an `int`. Size of the sliding window.
#' @param exclusion.zone a `numeric`. Size of the exclusion zone, based on window size (default is
#'   `1/2`). See details.
#' @param s.size a `numeric`. for anytime algorithm, represents the size (in observations) the
#'   random calculation will occur (default is `Inf`).
#' @param verbose an `int`. See details. (Default is `2`).
#'
#' @return Returns the matrix profile `mp` and profile index `pi`.
#' @export
#'
#' @family matrix profile computations
#' @seealso [mstomp()], [mstomp.par()]
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#' mp <- scrimp(toy_data$data[1:200,1], window.size = 30, verbose = 0)
#' \dontrun{
#' ref.data <- toy_data$data[,1]
#' query.data <- toy_data$data[,2]
#' # self similarity
#' mp <- scrimp(ref.data, window.size = 30, s.size = round(nrow(ref.data) * 0.1))
#' # join similarity
#' mp <- scrimp(ref.data, query.data, window.size = 30, s.size = round(nrow(query.data) * 0.1))
#' }

scrimp <- function(..., window.size, exclusion.zone = 1 / 2, s.size = Inf, verbose = 2) {
  args <- list(...)
  data <- args[[1]]
  if (length(args) > 1) {
    message("DISCLAIMER: This algorithm still in development by its authors.")
    message("Join similarity not implemented yet.")
    invisible(return(NULL))
    query <- args[[2]]
    exclusion.zone <- 0 # don't use exclusion zone for joins
  } else {
    query <- data
  }

  ## transform data into matrix
  if (is.vector(data)) {
    data <- as.matrix(data)
  }
  else if (is.matrix(data)) {
    if (ncol(data) > nrow(data)) {
      data <- t(data)
    }
  } else {
    stop("Error: Unknown type of data. Must be: a column matrix or a vector.", call. = FALSE)
  }

  if (is.vector(query)) {
    query <- as.matrix(query)
  } else if (is.matrix(query)) {
    if (ncol(query) > nrow(query)) {
      query <- t(query)
    }
  } else {
    stop("Error: Unknown type of query. Must be: a column matrix or a vector.", call. = FALSE)
  }

  exclusion.zone <- floor(window.size * exclusion.zone)
  data.size <- nrow(data)
  query.size <- nrow(query)
  matrix.profile.size <- data.size - window.size + 1
  num.queries <- query.size - window.size + 1

  ## check skip position
  skip.location <- rep(FALSE, matrix.profile.size)

  for (i in 1:matrix.profile.size) {
    if (any(is.na(data[i:(i + window.size - 1)])) || any(is.infinite(data[i:(i + window.size - 1)]))) {
      skip.location[i] <- TRUE
    }
  }

  data[is.na(data)] <- 0
  data[is.infinite(data)] <- 0

  query[is.na(query)] <- 0
  query[is.infinite(query)] <- 0

  if (query.size > data.size) {
    stop("Error: Query must be smaller or the same size as reference data.", call. = FALSE)
  }
  if (window.size > query.size / 2) {
    stop("Error: Time series is too short relative to desired window size.", call. = FALSE)
  }
  if (window.size < 4) {
    stop("Error: `window.size` must be at least 4.", call. = FALSE)
  }

  message("DISCLAIMER: This algorithm still in development by its authors.")

  matrix.profile <- matrix(Inf, matrix.profile.size, 1)
  left.matrix.profile <- right.matrix.profile <- matrix.profile
  profile.index <- matrix(-1, matrix.profile.size, 1)
  left.profile.index <- right.profile.index <- profile.index

  j <- 1
  order <- (exclusion.zone + 1):num.queries
  ssize <- min(s.size, length(order))
  order <- order[sample(1:length(order), size = ssize)]

  tictac <- Sys.time()

  if (verbose > 0) {
    pb <- utils::txtProgressBar(min = 1, max = ssize, style = 3, width = 80)
  }

  if (verbose > 0) {
    on.exit(close(pb))
  }
  if (verbose > 1) {
    on.exit(beep(sounds[[1]]), TRUE)
  }
  # anytime must return the result always
  on.exit(return(list(
    rmp = right.matrix.profile, rpi = right.profile.index,
    lmp = left.matrix.profile, lpi = left.profile.index,
    mp = Re(sqrt(as.complex(matrix.profile))), pi = profile.index
  )), TRUE)

  pre <- mass.pre(data, data.size, query, query.size, window.size = window.size)

  for (i in order) {
    j <- j + 1

    distance.profile <- diagonal.dist(
      data, i, data.size, window.size, num.queries, pre$data.mean,
      pre$data.sd
    )

    # distance.profile <- Re(sqrt(distance.profile))
    # distance.profile <- sqrt(abs(distance.profile))

    pos1 <- i:matrix.profile.size
    pos2 <- 1:(matrix.profile.size - i + 1)

    ind <- (matrix.profile[pos1] > distance.profile)
    profile.index[pos1[ind]] <- pos2[ind]
    matrix.profile[pos1[ind]] <- distance.profile[ind]
    ind <- (matrix.profile[pos2] > distance.profile)
    profile.index[pos2[ind]] <- pos1[ind]
    matrix.profile[pos2[ind]] <- distance.profile[ind]

    # matrix.profile[isSkip] <- Inf
    # profile.index[isSkip] <- 0

    # # anytime version
    # # left matrix.profile
    # ind <- (distance.profile[i:matrix.profile.size] < left.matrix.profile[i:matrix.profile.size])
    # ind <- c(rep(FALSE, (i - 1)), ind) # pad left
    # left.matrix.profile[ind] <- distance.profile[ind]
    # left.profile.index[which(ind)] <- i

    # # right matrix.profile
    # ind <- (distance.profile[1:i] < right.matrix.profile[1:i])
    # ind <- c(ind, rep(FALSE, matrix.profile.size - i)) # pad right
    # right.matrix.profile[ind] <- distance.profile[ind]
    # right.profile.index[which(ind)] <- i

    if (verbose > 0) {
      utils::setTxtProgressBar(pb, j)
    }
  }

  tictac <- Sys.time() - tictac

  if (verbose > 0) {
    message(sprintf("\nFinished in %.2f %s", tictac, units(tictac)))
  }

  # return() is at on.exit() function
}

#' SCRIMP NN algorithm
#'
#' @param data a `matrix` or a `vector`.
#' @param idx an `int`. The index of window to be computed.
#' @param data.size an `int`.
#' @param window.size an `int`.
#' @param mp.size an `int`.
#' @param data.mean result of `fast.movavg`.
#' @param data.sd result of `fast.movsd`.
#'
#' @return Returns the distance profile
#' @keywords internal
#' @noRd

diagonal.dist <- function(data, idx, data.size, window.size, mp.size, data.mean, data.sd) {
  data <- as.matrix(data)
  x.term <- matrix(1, mp.size - idx + 1, 1) *
    (t(data[idx:(idx + window.size - 1), 1, drop = FALSE]) %*% data[1:window.size, 1, drop = FALSE])[1]
  m.term <- data[idx:(mp.size - 1)] *
    data[1:(mp.size - idx)]
  a.term <- data[(idx + window.size):length(data)] *
    data[(window.size + 1):(data.size - idx + 1)]
  if (mp.size != idx) {
    x.term[2:length(x.term)] <- x.term[2:length(x.term)] - cumsum(m.term) + cumsum(a.term)
  }

  distance.profile <- (x.term - window.size * data.mean[idx:length(data.mean)] * data.mean[1:(mp.size - idx + 1)]) /
    (window.size * data.sd[idx:length(data.sd)] * data.sd[1:(mp.size - idx + 1)])
  distance.profile <- 2 * window.size * (1 - distance.profile)

  return(distance.profile)
}
