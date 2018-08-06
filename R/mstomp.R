#' Multivariate STOMP algorithm
#'
#' No description
#'
#' No details
#'
#' @param data a matrix, where each colums is a time series (dimention). Can accept Lists and data.frames too.
#' @param window.size size of the sliding window
#' @param must.dim which dimentions to forcibly include (default is NULL)
#' @param exc.dim which dimentions to exclude (default is NULL)
#'
#' @return The matrix profile and profile index
#' @export
#'
#' @seealso [stamp()]
#' @references Yeh CM, Kavantzas N, Keogh E. Matrix Profile VI : Meaningful Multidimensional Motif Discovery.
#' [https://sites.google.com/view/mstamp/]
#'
#' @examples
#' \dontrun{
#' mp <- mstomp(data, 30)
#' mp <- mstomp(data, 30, must.dim = c(1, 2))
#' mp <- mstomp(data, 30, exc.dim = c(2,3))
#' }
mstomp <- function(data, window.size, must.dim = NULL, exc.dim = NULL, exclusion.zone = 1 / 2) {

  eps <- .Machine$double.eps^0.5

  ## get various length
  exclusion.zone <- floor(window.size * exclusion.zone)

  ## transform data list into matrix
  if (is.list(data)) {
    data.size <- length(data[[1]])
    n.dim <- length(data)

    for (i in 1:n.dim) {
      len <- length(data[[i]])
      # Fix TS size with NaN
      if (len < data.size) {
        data[[i]] <- c(data[[i]], rep(NA, data.size - len))
      }
    }
    # transform data into matrix (each column is a TS)
    data <- sapply(data, cbind)
  } else if (is.matrix(data) || is.data.frame(data)) {
    if (is.data.frame(data)) {
      data <- as.matrix(data)
    } # just to be uniform
    if (ncol(data) > nrow(data)) {
      data <- t(data)
    }
    data.size <- nrow(data)
    n.dim <- ncol(data)
  } else if (is.vector(data)) {
    data.size <- length(data)
    n.dim <- 1
    # transform data into 1-col matrix
    data <- as.matrix(data) # just to be uniform
  } else {
    stop("Unknown type of data. Must be: matrix, data.frame, vector or list")
  }

  matrix.profile.size <- data.size - window.size + 1

  ## check input
  if (window.size > data.size / 2) {
    stop("Error: Time series is too short relative to desired subsequence length")
  }
  if (window.size < 4) {
    stop("Error: Subsequence length must be at least 4")
  }
  if (any(must.dim > n.dim)) {
    stop("Error: The must have dimension must be less then the total dimension")
  }
  if (any(exc.dim > n.dim)) {
    stop("Error: The exclusion dimension must be less then the total dimension")
  }
  if (length(intersect(must.dim, exc.dim)) > 0) {
    stop("Error: The same dimension is presented in both the exclusion dimension and must have dimension")
  }

  ## check skip position
  n.exc <- length(exc.dim)
  n.must <- length(must.dim)
  mask.exc <- rep(FALSE, n.dim)
  mask.exc[exc.dim] <- TRUE
  skip.location <- rep(FALSE, matrix.profile.size)

  for (i in 1:matrix.profile.size) {
    if (any(is.na(data[i:(i + window.size - 1), !mask.exc])) || any(is.infinite(data[i:(i + window.size - 1), !mask.exc]))) {
      skip.location[i] <- TRUE
    }
  }

  data[is.na(data)] <- 0
  data[is.infinite(data)] <- 0

  pb <- txtProgressBar(min = 0, max = matrix.profile.size, style = 3)
  on.exit(close(pb))
  on.exit(beep(), TRUE)

  ## initialization
  data.fft <- matrix(0, (window.size + data.size), n.dim)
  data.mean <- matrix(0, matrix.profile.size, n.dim)
  data.sd <- matrix(0, matrix.profile.size, n.dim)
  first.profile <- matrix(0, matrix.profile.size, n.dim)

  for (i in 1:n.dim) {
    nnPre <- mass.pre(data[, i], data.size, window.size = window.size)
    data.fft[, i] <- nnPre$data.fft
    data.mean[, i] <- nnPre$data.mean
    data.sd[, i] <- nnPre$data.sd
    mstomp <- mass(data.fft[, i], data[1:window.size, i], data.size, window.size, data.mean[, i], data.sd[, i], data.mean[1, i], data.sd[1, i])
    first.profile[, i] <- mstomp$last.product
  }

  ## compute the matrix profile
  matrix.profile <- matrix(0, matrix.profile.size, n.dim)
  profile.index <- matrix(0, matrix.profile.size, n.dim)
  distance.profile <- matrix(0, matrix.profile.size, n.dim)
  last.product <- matrix(0, matrix.profile.size, n.dim)
  drop.value <- matrix(0, 1, n.dim)
  for (i in 1:matrix.profile.size) {
    # compute the distance profile
    setTxtProgressBar(pb, i)

    query <- as.matrix(data[i:(i + window.size - 1), ])

    if (i == 1) {
      for (j in 1:n.dim) {
        mstomp <- mass(data.fft[, j], query[, j], data.size, window.size, data.mean[, j], data.sd[, j], data.mean[i, j], data.sd[i, j])
        distance.profile[, j] <- mstomp$distance.profile
        last.product[, j] <- mstomp$last.product
      }
    } else {
      rep.drop.val <- kronecker(matrix(1, matrix.profile.size - 1, 1), t(drop.value))
      rep.query <- kronecker(matrix(1, matrix.profile.size - 1, 1), t(query[window.size, ]))

      last.product[2:(data.size - window.size + 1), ] <- last.product[1:(data.size - window.size), ] -
        data[1:(data.size - window.size), ] * rep.drop.val +
        data[(window.size + 1):data.size, ] * rep.query


      last.product[1, ] <- first.profile[i, ]

      distance.profile <- 2 * (window.size - (last.product - window.size * data.mean * kronecker(matrix(1, matrix.profile.size, 1), t(data.mean[i, ]))) /
        (data.sd * kronecker(matrix(1, matrix.profile.size, 1), t(data.sd[i, ]))))
    }

    distance.profile <- Re(distance.profile)
    drop.value <- query[1, ]

    # apply exclusion zone
    exc.st <- max(1, i - exclusion.zone)
    exc.ed <- min(matrix.profile.size, i + exclusion.zone)
    distance.profile[exc.st:exc.ed, ] <- Inf
    distance.profile[data.sd < eps] <- Inf
    if (skip.location[i] || any(data.sd[i, !mask.exc] < eps)) {
      distance.profile[] <- Inf
    }
    distance.profile[skip.location, ] <- Inf

    # apply dimension "must have" and "exclusion"
    distance.profile[, exc.dim] <- Inf

    if (n.must > 0) {
      mask.must <- rep(FALSE, n.must)
      mask.must[must.dim] <- TRUE
      dist.pro.must <- distance.profile[, mask.must]
      distance.profile[, mask.must] <- -Inf
    }

    if (n.dim > 1) {
      dist.pro.sort <- t(apply(distance.profile, 1, sort))
    } # sort by row, put all -Inf to the first columns
    else {
      dist.pro.sort <- distance.profile
    }

    if (n.must > 0) {
      dist.pro.sort[, 1:n.must] <- dist.pro.must
    }

    # figure out and store the nearest neighbor
    dist.pro.cum <- rep(0, matrix.profile.size)
    dist.pro.merg <- rep(0, matrix.profile.size)
    for (j in (max(1, n.must):(n.dim - n.exc))) {
      dist.pro.cum <- dist.pro.cum + dist.pro.sort[, j]
      dist.pro.merg[] <- dist.pro.cum / j
      min.idx <- which.min(dist.pro.merg)
      min.val <- dist.pro.merg[min.idx]
      matrix.profile[i, j] <- min.val
      profile.index[i, j] <- min.idx
    }
  }

  ## remove bad k setting in the returned matrix
  matrix.profile <- sqrt(matrix.profile)
  if (n.must > 0) {
    matrix.profile[, 1:(n.must - 1)] <- NA
  }
  if (n.exc > 0) {
    matrix.profile[, (n.dim - n.exc + 1):length(matrix.profile)] <- NA
  }
  if (n.must > 0) {
    profile.index[, 1:(n.must - 1)] <- NA
  }
  if (n.exc > 0) {
    profile.index[, (n.dim - n.exc + 1):length(profile.index)] <- NA
  }

  return(list(mp = matrix.profile, pi = profile.index))
}
