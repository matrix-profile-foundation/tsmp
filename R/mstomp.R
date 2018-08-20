#' Multivariate STOMP algorithm
#'
#' Computes the Matrix Profile and Profile Index for Multivariate Time Series.
#'
#' The Matrix Profile, has the potential to revolutionize time series data mining because of its generality, versatility, simplicity and scalability. In particular it has implications for time series motif discovery, time series joins, shapelet discovery (classification), density estimation, semantic segmentation, visualization, rule discovery, clustering etc.
#' The MSTOMP computes the Matrix Profile and Profile Index for Multivariate Time Series that is meaningful for multidimensional MOTIF discovery. It uses the STOMP algorithm that is faster than STAMP but lacks its anytime property.
#'
#' Although this functions handles Multivariate Time Series, it can also be used to handle Univariate Time Series.
#' `verbose` changes how much information is printed by this function; `0` means nothing, `1` means text, `2` means text and sound.
#'
#' @param data a `matrix` of `numeric`, where each column is a time series. Accepts `vector` (see details), `list` and `data.frame` too.
#' @param window.size an `int` with the size of the sliding window.
#' @param must.dim an `int` or `vector` of which dimensions to forcibly include (default is `NULL`).
#' @param exc.dim an `int` or `vector` of which dimensions to exclude (default is `NULL`).
#' @param exclusion.zone a `numeric`. Size of the exclusion zone, based on query size (default is `1/2`).
#' @param verbose an `int`. See details. (Default is `2`).
#'
#' @return Returns the matrix profile `mp` and profile index `pi`.
#' It also returns the left and right matrix profile `lmp`, `rmp` and profile index `lpi`, `rpi` that may be used to detect Time Series Chains (Yan Zhu 2018).
#' @export
#'
#' @family mstomp
#' @seealso [stamp()], [stamp.par()], [mstomp.par()]
#' @references 1. Yeh CM, Kavantzas N, Keogh E. Matrix Profile VI : Meaningful Multidimensional Motif Discovery.
#' @references 2. Zhu Y, Imamura M, Nikovski D, Keogh E. Matrix Profile VII: Time Series Chains: A New Primitive for Time Series Data Mining. Knowl Inf Syst. 2018 Jun 2;1–27.
#' @references Website: <https://sites.google.com/view/mstamp/>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#' # using all dimensions
#' mp <- mstomp(toy_data$data[1:200,], 30, verbose = 0)
#' \dontrun{
#' # force using dimensions 1 and 2
#' mp <- mstomp(toy_data$data[1:200,], 30, must.dim = c(1, 2))
#' # exclude dimensions 2 and 3
#' mp <- mstomp(toy_data$data[1:200,], 30, exc.dim = c(2, 3))
#' }

mstomp <- function(data, window.size, must.dim = NULL, exc.dim = NULL, exclusion.zone = 1 / 2, verbose = 2) {
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

  if (verbose > 0) {
    pb <- utils::txtProgressBar(min = 0, max = matrix.profile.size, style = 3, width = 80)
    on.exit(close(pb))
  }
  if (verbose > 1) {
    on.exit(beep(sounds[[1]]), TRUE)
  }

  ## initialization
  data.fft <- matrix(0, (window.size + data.size), n.dim)
  data.mean <- matrix(0, matrix.profile.size, n.dim)
  data.sd <- matrix(0, matrix.profile.size, n.dim)
  first.product <- matrix(0, matrix.profile.size, n.dim)

  for (i in 1:n.dim) {
    nnPre <- mass.pre(data[, i], data.size, window.size = window.size)
    data.fft[, i] <- nnPre$data.fft
    data.mean[, i] <- nnPre$data.mean
    data.sd[, i] <- nnPre$data.sd
    mstomp <- mass(data.fft[, i], data[1:window.size, i], data.size, window.size, data.mean[, i], data.sd[, i], data.mean[1, i], data.sd[1, i])
    first.product[, i] <- mstomp$last.product
  }

  tictac <- Sys.time()
  ## compute the matrix profile
  matrix.profile <- matrix(0, matrix.profile.size, n.dim)
  profile.index <- matrix(0, matrix.profile.size, n.dim)
  left.matrix.profile <- matrix(Inf, matrix.profile.size, n.dim)
  left.profile.index <- matrix(-1, matrix.profile.size, n.dim)
  right.matrix.profile <- matrix(Inf, matrix.profile.size, n.dim)
  right.profile.index <- matrix(-1, matrix.profile.size, n.dim)
  distance.profile <- matrix(0, matrix.profile.size, n.dim)
  last.product <- matrix(0, matrix.profile.size, n.dim)
  drop.value <- matrix(0, 1, n.dim)

  for (i in 1:matrix.profile.size) {
    # compute the distance profile
    if (verbose > 0) {
      utils::setTxtProgressBar(pb, i)
    }

    query <- as.matrix(data[i:(i + window.size - 1), ])

    if (i == 1) {
      for (j in 1:n.dim) {
        mstomp <- mass(data.fft[, j], query[, j], data.size, window.size, data.mean[, j], data.sd[, j], data.mean[i, j], data.sd[i, j])
        distance.profile[, j] <- mstomp$distance.profile
        last.product[, j] <- mstomp$last.product
      }
    } else {
      rep.drop.value <- kronecker(matrix(1, matrix.profile.size - 1, 1), t(drop.value))
      rep.query <- kronecker(matrix(1, matrix.profile.size - 1, 1), t(query[window.size, ]))

      last.product[2:(data.size - window.size + 1), ] <- last.product[1:(data.size - window.size), ] -
        data[1:(data.size - window.size), ] * rep.drop.value +
        data[(window.size + 1):data.size, ] * rep.query


      last.product[1, ] <- first.product[i, ]

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
      mask.must <- rep(FALSE, n.dim)
      mask.must[must.dim] <- TRUE
      dist.pro.must <- distance.profile[, mask.must]
      distance.profile[, mask.must] <- -Inf
    }

    if (n.dim > 1) {
      dist.pro.sort <- t(apply(distance.profile, 1, sort))
    } # sort by row, put all -Inf to the first column
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

      # left matrix.profile
      if (i > (exclusion.zone + 1)) {
        min.idx <- which.min(dist.pro.merg[1:(i - exclusion.zone)])
        min.val <- dist.pro.merg[min.idx]
        left.matrix.profile[i, j] <- min.val
        left.profile.index[i, j] <- min.idx
      }

      # right matrix.profile
      if (i < (matrix.profile.size - exclusion.zone)) {
        min.idx <- which.min(dist.pro.merg[(i + exclusion.zone):matrix.profile.size]) + i + exclusion.zone - 1
        min.val <- dist.pro.merg[min.idx]
        right.matrix.profile[i, j] <- min.val
        right.profile.index[i, j] <- min.idx
      }

      # normal matrix.profile
      min.idx <- which.min(dist.pro.merg)
      min.val <- dist.pro.merg[min.idx]
      matrix.profile[i, j] <- min.val
      profile.index[i, j] <- min.idx
    }
  }

  matrix.profile <- sqrt(matrix.profile)
  right.matrix.profile <- sqrt(right.matrix.profile)
  left.matrix.profile <- sqrt(left.matrix.profile)

  ## remove bad k setting in the returned matrix
  if (n.must > 1) {
    matrix.profile[, 1:(n.must - 1)] <- NA
    right.matrix.profile[, 1:(n.must - 1)] <- NA
    left.matrix.profile[, 1:(n.must - 1)] <- NA
  }
  if (n.exc > 0) {
    matrix.profile[, (n.dim - n.exc + 1):n.dim] <- NA
    right.matrix.profile[, (n.dim - n.exc + 1):n.dim] <- NA
    left.matrix.profile[, (n.dim - n.exc + 1):n.dim] <- NA
  }
  if (n.must > 1) {
    profile.index[, 1:(n.must - 1)] <- NA
    right.profile.index[, 1:(n.must - 1)] <- NA
    left.profile.index[, 1:(n.must - 1)] <- NA
  }
  if (n.exc > 0) {
    profile.index[, (n.dim - n.exc + 1):n.dim] <- NA
    right.profile.index[, (n.dim - n.exc + 1):n.dim] <- NA
    left.profile.index[, (n.dim - n.exc + 1):n.dim] <- NA
  }

  tictac <- Sys.time() - tictac

  if (verbose > 0) {
    message(sprintf("\nFinished in %.2f %s", tictac, units(tictac)))
  }

  return(list(
    rmp = right.matrix.profile, rpi = right.profile.index,
    lmp = left.matrix.profile, lpi = left.profile.index,
    mp = matrix.profile, pi = profile.index
  ))
}
