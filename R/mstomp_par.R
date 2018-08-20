#' Multivariate STOMP algorithm Parallel version
#'
#' Computes the Matrix Profile and Profile Index for Multivariate Time Series.
#'
#' The Matrix Profile, has the potential to revolutionize time series data mining because of its generality, versatility, simplicity and scalability. In particular it has implications for time series motif discovery, time series joins, shapelet discovery (classification), density estimation, semantic segmentation, visualization, rule discovery, clustering etc.
#' The MSTOMP computes the Matrix Profile and Profile Index for Multivariate Time Series that is meaningful for multidimensional MOTIF discovery. It uses the STOMP algorithm that is faster than STAMP but lacks its anytime property.
#'
#' Although this functions handles Multivariate Time Series, it can also be used to handle Univariate Time Series.
#' `verbose` changes how much information is printed by this function; `0` means nothing, `1` means text, `2` means text and sound.
#'
#' @param data a `matrix` of `numeric`, where each colums is a time series. Accepts `vector` (see details), `list` and `data.frame` too.
#' @param window.size an `int`. Size of the sliding window.
#' @param must.dim an `int` or `vector` of which dimensions to forcibly include (default is `NULL`).
#' @param exc.dim an `int` or `vector` of which dimensions to exclude (default is `NULL`).
#' @param exclusion.zone a `numeric`. Size of the exclusion zone, based on query size (default is `1/2`). See details.
#' @param verbose an `int`. See details. (Default is `2`).
#' @param n.workers an `int`. Number of workers for parallel. (Default is `2`).
#'
#' @return Returns the matrix profile `mp` and profile index `pi`.
#' It also returns the left and right matrix profile `lmp`, `rmp` and profile index `lpi`, `rpi` that may be used to detect Time Series Chains (Yan Zhu 2018).
#' @export
#'
#' @family mstomp
#' @seealso [stamp()], [stamp.par()], [mstomp()]
#' @references 1. Yeh CM, Kavantzas N, Keogh E. Matrix Profile VI : Meaningful Multidimensional Motif Discovery.
#' @references 2. Zhu Y, Imamura M, Nikovski D, Keogh E. Matrix Profile VII: Time Series Chains: A New Primitive for Time Series Data Mining. Knowl Inf Syst. 2018 Jun 2;1–27.
#' @references Website: <https://sites.google.com/view/mstamp/>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#' # using all dimensions
#' Sys.sleep(1) # sometimes sleep is needed if you run parallel multiple times in a row
#' mp <- mstomp.par(toy_data$data[1:100,], 30, verbose = 0)
#' @import beepr doSNOW foreach parallel

mstomp.par <- function(data, window.size, must.dim = NULL, exc.dim = NULL, exclusion.zone = 1 / 2, verbose = 2, n.workers = 2) {
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

  cores <- min(max(2, n.workers), parallel::detectCores())

  # SNOW package
  if (verbose > 0) {
    progress <- function(n) utils::setTxtProgressBar(pb, n)
  }
  else {
    progress <- function(n) return(invisible(TRUE))
  }
  opts <- list(progress = progress)

  cl <- parallel::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  on.exit(parallel::stopCluster(cl))
  if (verbose > 0) {
    on.exit(close(pb), TRUE)
  }
  if (verbose > 1) {
    on.exit(beepr::beep(10), TRUE)
  }

  ## initialize variable
  per.work <- max(10, ceiling(matrix.profile.size / 100))
  n.work <- floor(matrix.profile.size / per.work)
  idx.work <- list()

  for (i in 1:n.work) {
    idx.st <- (i - 1) * per.work + 1
    if (i == n.work) {
      idx.ed <- matrix.profile.size
    } else {
      idx.ed <- i * per.work
    }
    idx.work[[i]] <- idx.st:idx.ed
  }

  tictac <- Sys.time()

  if (verbose > 0) {
    pb <- utils::txtProgressBar(min = 0, max = n.work, style = 3, width = 80)
  }

  i <- NULL # CRAN NOTE fix
  `%dopar%` <- foreach::`%dopar%` # CRAN NOTE fix

  ## compute the matrix profile
  batch <- foreach(
    i = 1:n.work,
    .verbose = FALSE,
    .inorder = FALSE,
    .multicombine = TRUE,
    .options.snow = opts,
    # .combine = combiner,
    # .errorhandling = 'remove',
    .export = c("mass")
  ) %dopar% {
    pro.muls <- matrix(0, length(idx.work[[i]]), n.dim)
    pro.idxs <- matrix(0, length(idx.work[[i]]), n.dim)
    pro.muls.right <- matrix(Inf, length(idx.work[[i]]), n.dim)
    pro.idxs.right <- matrix(-1, length(idx.work[[i]]), n.dim)
    pro.muls.left <- matrix(Inf, length(idx.work[[i]]), n.dim)
    pro.idxs.left <- matrix(-1, length(idx.work[[i]]), n.dim)
    dist.pro <- matrix(0, matrix.profile.size, n.dim)
    last.product <- matrix(0, matrix.profile.size, n.dim)
    drop.value <- matrix(0, 1, n.dim)

    for (j in 1:length(idx.work[[i]])) {
      idx <- idx.work[[i]][j]

      query <- as.matrix(data[idx:(idx + window.size - 1), ])

      if (j == 1) {
        for (k in 1:n.dim) {
          mstomp <- mass(data.fft[, k], query[, k], data.size, window.size, data.mean[, k], data.sd[, k], data.mean[idx, k], data.sd[idx, k])
          dist.pro[, k] <- mstomp$distance.profile
          last.product[, k] <- mstomp$last.product
        }
      } else {
        rep.drop.value <- kronecker(matrix(1, matrix.profile.size - 1, 1), t(drop.value))
        rep.query <- kronecker(matrix(1, matrix.profile.size - 1, 1), t(query[window.size, ]))

        last.product[2:(data.size - window.size + 1), ] <- last.product[1:(data.size - window.size), ] -
          data[1:(data.size - window.size), ] * rep.drop.value +
          data[(window.size + 1):data.size, ] * rep.query


        last.product[1, ] <- first.product[idx, ]

        dist.pro <- 2 * (window.size - (last.product - window.size * data.mean * kronecker(matrix(1, matrix.profile.size, 1), t(data.mean[idx, ]))) /
          (data.sd * kronecker(matrix(1, matrix.profile.size, 1), t(data.sd[idx, ]))))
      }

      dist.pro <- Re(dist.pro)
      # dist.pro <- max(dist.pro, 0)
      drop.value <- query[1, ]

      # apply exclusion zone
      exc.zone.st <- max(1, idx - exclusion.zone)
      exc.zone.ed <- min(matrix.profile.size, idx + exclusion.zone)
      dist.pro[exc.zone.st:exc.zone.ed, ] <- Inf
      dist.pro[data.sd < eps] <- Inf
      if (skip.location[idx] || any(data.sd[idx, !mask.exc] < eps)) {
        dist.pro[] <- Inf
      }
      dist.pro[skip.location, ] <- Inf

      # apply dimension "must have" and "exclusion"
      dist.pro[, exc.dim] <- Inf

      if (n.must > 0) {
        mask.must <- rep(FALSE, n.dim)
        mask.must[must.dim] <- TRUE
        dist.pro.must <- dist.pro[, mask.must]
        dist.pro[, mask.must] <- -Inf
      }

      # figure out and store the nearest neighbor
      if (n.dim > 1) {
        dist.pro.sort <- t(apply(dist.pro, 1, sort))
      } # sort by row, left to right
      else {
        dist.pro.sort <- dist.pro
      }

      if (n.must > 0) {
        dist.pro.sort[, 1:n.must] <- dist.pro.must
      }
      dist.pro.cum <- rep(0, matrix.profile.size)
      dist.pro.merg <- rep(0, matrix.profile.size)

      for (k in (max(1, n.must):(n.dim - n.exc))) {
        dist.pro.cum <- dist.pro.cum + dist.pro.sort[, k]
        dist.pro.merg[] <- dist.pro.cum / k
        # left matrix.profile
        if (idx > (exclusion.zone + 1)) {
          min.idx <- which.min(dist.pro.merg[1:(idx - exclusion.zone)])
          min.val <- dist.pro.merg[min.idx]
          pro.muls.left[j, k] <- min.val
          pro.idxs.left[j, k] <- min.idx
        }

        # right matrix.profile
        if (idx < (matrix.profile.size - exclusion.zone)) {
          min.idx <- which.min(dist.pro.merg[(idx + exclusion.zone):matrix.profile.size]) + idx + exclusion.zone - 1
          min.val <- dist.pro.merg[min.idx]
          pro.muls.right[j, k] <- min.val
          pro.idxs.right[j, k] <- min.idx
        }

        # normal matrix.profile
        min.idx <- which.min(dist.pro.merg)
        min.val <- dist.pro.merg[min.idx]
        pro.muls[j, k] <- min.val
        pro.idxs[j, k] <- min.idx
      }
    }

    pro.muls <- sqrt(pro.muls)
    pro.muls.left <- sqrt(pro.muls.left)
    pro.muls.right <- sqrt(pro.muls.right)

    res <- list(pro.muls = pro.muls, pro.idxs = pro.idxs, pro.muls.left = pro.muls.left, pro.idxs.left = pro.idxs.left, pro.muls.right = pro.muls.right, pro.idxs.right = pro.idxs.right, idx = i)

    res
  }

  matrix.profile <- matrix(0, matrix.profile.size, n.dim)
  profile.index <- matrix(0, matrix.profile.size, n.dim)
  left.matrix.profile <- matrix(Inf, matrix.profile.size, n.dim)
  left.profile.index <- matrix(-1, matrix.profile.size, n.dim)
  right.matrix.profile <- matrix(Inf, matrix.profile.size, n.dim)
  right.profile.index <- matrix(-1, matrix.profile.size, n.dim)

  for (i in 1:length(batch)) {
    left.profile.index[idx.work[[batch[[i]]$idx]], ] <- batch[[i]]$pro.idxs.left
    left.matrix.profile[idx.work[[batch[[i]]$idx]], ] <- batch[[i]]$pro.muls.left
    right.profile.index[idx.work[[batch[[i]]$idx]], ] <- batch[[i]]$pro.idxs.right
    right.matrix.profile[idx.work[[batch[[i]]$idx]], ] <- batch[[i]]$pro.muls.right
    profile.index[idx.work[[batch[[i]]$idx]], ] <- batch[[i]]$pro.idxs
    matrix.profile[idx.work[[batch[[i]]$idx]], ] <- batch[[i]]$pro.muls
  }

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
