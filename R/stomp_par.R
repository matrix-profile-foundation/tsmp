#' Univariate STOMP algorithm
#'
#' Computes the Matrix Profile and Profile Index for Univariate Time Series.
#'
#' @details
#' The Matrix Profile, has the potential to revolutionize time series data mining because of its
#' generality, versatility, simplicity and scalability. In particular it has implications for time
#' series motif discovery, time series joins, shapelet discovery (classification), density
#' estimation, semantic segmentation, visualization, rule discovery, clustering etc. `verbose`
#' changes how much information is printed by this function; `0` means nothing, `1` means text, `2`
#' means text and sound. `exclusion.zone` is used to avoid trivial matches; if a query data is
#' provided (join similarity), this parameter is ignored.
#'
#' @param ... a `matrix` or a `vector`. If a second time series is supplied it will be a join matrix
#'   profile.
#' @param window.size an `int`. Size of the sliding window.
#' @param exclusion.zone a `numeric`. Size of the exclusion zone, based on window size (default is
#'   `1/2`). See details.
#' @param verbose an `int`. See details. (Default is `2`).
#' @param n.workers an `int`. Number of workers for parallel. (Default is `2`).
#'
#' @return Returns the matrix profile `mp` and profile index `pi`. It also returns the left and
#'   right matrix profile `lmp`, `rmp` and profile index `lpi`, `rpi` that may be used to detect
#'   Time Series Chains (Yan Zhu 2018).
#' @export
#'
#' @family matrix profile computations
#' @seealso [stamp()], [stamp.par()]; [mstomp()], [mstomp.par()] for multivariate analysis.
#' @references * Zhu Y, Zimmerman Z, Senobari NS, Yeh CM, Funning G. Matrix Profile II : Exploiting
#'   a Novel Algorithm and GPUs to Break the One Hundred Million Barrier for Time Series Motifs and
#'   Joins. Icdm. 2016 Jan 22;54(1):739–48.
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#' mp <- stomp.par(toy_data$data[1:200,1], window.size = 30, verbose = 0)
#' \dontrun{
#' ref.data <- toy_data$data[,1]
#' query.data <- toy_data$data[,2]
#' # self similarity
#' mp <- stomp.par(ref.data, window.size = 30)
#' # join similarity
#' mp <- stomp.par(ref.data, query.data, window.size = 30)
#' }
stomp.par <- function(..., window.size, exclusion.zone = 1 / 2, verbose = 2, n.workers = 2) {
  args <- list(...)
  data <- args[[1]]
  if (length(args) > 1) {
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

  if (query.size > data.size) {
    stop("Error: Query must be smaller or the same size as reference data.", call. = FALSE)
  }
  if (window.size > query.size / 2) {
    stop("Error: Time series is too short relative to desired window size.", call. = FALSE)
  }
  if (window.size < 4) {
    stop("Error: `window.size` must be at least 4.", call. = FALSE)
  }

  data.fft <- matrix(0, (window.size + data.size), 1)
  data.mean <- matrix(0, matrix.profile.size, 1)
  data.sd <- matrix(0, matrix.profile.size, 1)
  first.product <- matrix(0, num.queries, 1)

  # forward
  nnpre <- mass.pre(data, data.size, query, query.size, window.size = window.size)
  data.fft <- nnpre$data.fft
  data.mean <- nnpre$data.mean
  data.sd <- nnpre$data.sd
  query.mean <- nnpre$query.mean
  query.sd <- nnpre$query.sd

  # reverse
  # This is needed to handle with the join similarity.
  rnnpre <- mass.pre(query, query.size, data, data.size, window.size = window.size)
  rdata.fft <- rnnpre$data.fft
  rdata.mean <- rnnpre$data.mean
  rdata.sd <- rnnpre$data.sd
  rquery.mean <- rnnpre$query.mean
  rquery.sd <- rnnpre$query.sd
  rnn <- mass(
    rdata.fft, data[1:window.size], query.size, window.size, rdata.mean,
    rdata.sd, rquery.mean[1], rquery.sd[1]
  )

  first.product[, 1] <- rnn$last.product

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
    on.exit(beep(sounds[[1]]), TRUE)
  }

  ## seperate index into different job
  per.work <- max(10, ceiling(num.queries / 100))
  n.work <- floor(num.queries / per.work)
  idx.work <- list()

  for (i in 1:n.work) {
    idx.st <- (i - 1) * per.work + 1
    if (i == n.work) {
      idx.ed <- num.queries
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
    .export = "mass"
  ) %dopar% {
    work.len <- length(idx.work[[i]])
    pro.muls <- matrix(Inf, matrix.profile.size, 1)
    pro.idxs <- matrix(-1, matrix.profile.size, 1)
    pro.muls.right <- matrix(Inf, matrix.profile.size, 1)
    pro.idxs.right <- matrix(-1, matrix.profile.size, 1)
    pro.muls.left <- matrix(Inf, matrix.profile.size, 1)
    pro.idxs.left <- matrix(-1, matrix.profile.size, 1)
    dist.pro <- matrix(0, matrix.profile.size, 1)
    last.product <- matrix(0, matrix.profile.size, 1)
    drop.value <- matrix(0, 1, 1)

    for (j in 1:work.len) {
      # compute the distance profile
      idx.st <- idx.work[[i]][1]
      idx.ed <- idx.work[[i]][work.len]
      idx <- idx.work[[i]][j]

      query.window <- as.matrix(query[idx:(idx + window.size - 1), 1])

      if (j == 1) {
        nn <- mass(
          data.fft, query.window, data.size, window.size, data.mean, data.sd,
          query.mean[idx], query.sd[idx]
        )
        dist.pro[, 1] <- nn$distance.profile
        last.product[, 1] <- nn$last.product
      } else {
        last.product[2:(data.size - window.size + 1), 1] <- last.product[1:(data.size - window.size), 1] -
          data[1:(data.size - window.size), 1] * drop.value +
          data[(window.size + 1):data.size, 1] * query.window[window.size, 1]
        last.product[1, 1] <- first.product[idx, 1]
        dist.pro <- 2 * (window.size - (last.product - window.size * data.mean * query.mean[idx]) /
          (data.sd * query.sd[idx]))
      }

      dist.pro <- Re(sqrt(dist.pro))
      drop.value <- query.window[1, 1]

      # apply exclusion zone
      if (exclusion.zone > 0) {
        exc.st <- max(1, idx - exclusion.zone)
        exc.ed <- min(matrix.profile.size, idx + exclusion.zone)
        dist.pro[exc.st:exc.ed, 1] <- Inf
        dist.pro[data.sd < vars()$eps] <- Inf
      }

      # left matrix.profile
      ind <- (dist.pro[idx:matrix.profile.size] < pro.muls.left[idx:matrix.profile.size])
      ind <- c(rep(FALSE, (idx - 1)), ind) # pad left
      pro.muls.left[ind] <- dist.pro[ind]
      pro.idxs.left[which(ind)] <- idx

      # right matrix.profile
      ind <- (dist.pro[1:idx] < pro.muls.right[1:idx])
      ind <- c(ind, rep(FALSE, matrix.profile.size - idx)) # pad right
      pro.muls.right[ind] <- dist.pro[ind]
      pro.idxs.right[which(ind)] <- idx

      # normal matrix.profile
      ind <- (dist.pro < pro.muls)
      pro.muls[ind] <- dist.pro[ind]
      pro.idxs[which(ind)] <- idx
    }

    res <- list(
      pro.muls = pro.muls, pro.idxs = pro.idxs,
      pro.muls.left = pro.muls.left, pro.idxs.left = pro.idxs.left,
      pro.muls.right = pro.muls.right, pro.idxs.right = pro.idxs.right
    )

    res
  }

  matrix.profile <- matrix(Inf, matrix.profile.size, 1)
  profile.index <- matrix(-1, matrix.profile.size, 1)
  left.matrix.profile <- matrix(Inf, matrix.profile.size, 1)
  left.profile.index <- matrix(-1, matrix.profile.size, 1)
  right.matrix.profile <- matrix(Inf, matrix.profile.size, 1)
  right.profile.index <- matrix(-1, matrix.profile.size, 1)

  for (i in 1:length(batch)) {
    ind <- (batch[[i]]$pro.muls < matrix.profile)
    matrix.profile[ind] <- batch[[i]]$pro.muls[ind]
    profile.index[ind] <- batch[[i]]$pro.idxs[ind]

    ind <- (batch[[i]]$pro.muls.left < left.matrix.profile)
    left.matrix.profile[ind] <- batch[[i]]$pro.muls.left[ind]
    left.profile.index[ind] <- batch[[i]]$pro.idxs.left[ind]

    ind <- (batch[[i]]$pro.muls.right < right.matrix.profile)
    right.matrix.profile[ind] <- batch[[i]]$pro.muls.right[ind]
    right.profile.index[ind] <- batch[[i]]$pro.idxs.right[ind]
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
