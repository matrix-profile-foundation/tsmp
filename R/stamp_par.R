#' Anytime univariate STAMP algorithm Parallel version
#'
#' Computes the best so far Matrix Profile and Profile Index for Univariate Time Series.
#'
#' The Matrix Profile, has the potential to revolutionize time series data mining because of its
#' generality, versatility, simplicity and scalability. In particular it has implications for time
#' series motif discovery, time series joins, shapelet discovery (classification), density
#' estimation, semantic segmentation, visualization, rule discovery, clustering etc. The anytime
#' STAMP computes the Matrix Profile and Profile Index in such manner that it can be stopped before
#' its complete calculation and return the best so far results allowing ultra-fast approximate
#' solutions. `verbose` changes how much information is printed by this function; `0` means nothing,
#' `1` means text, `2` means text and sound. `exclusion.zone` is used to avoid  trivial matches; if
#' a query data is provided (join similarity), this parameter is ignored.
#'
#' @param ... a `matrix` or a `vector`. If a second time series is supplied it will be a join matrix
#'   profile.
#' @param window.size an `int`. Size of the sliding window.
#' @param exclusion.zone a `numeric`. Size of the exclusion zone, based on query size (default is
#'   `1/2`). See details.
#' @param s.size a `numeric`. for anytime algorithm, represents the size (in observations) the
#'   random calculation will occur (default is `Inf`).
#' @param n.workers an `int`. Number of workers for parallel. (Default is `2`).
#' @param verbose an `int`. See details. (Default is `2`).
#'
#' @return Returns the matrix profile `mp` and profile index `pi`. It also returns the left and
#'   right matrix profile `lmp`, `rmp` and profile index `lpi`, `rpi` that may be used to detect
#'   Time Series Chains (Yan Zhu 2018).
#' @export
#'
#' @family Stamp
#' @seealso [mstomp()], [mstomp.par()]
#' @references * Yeh CCM, Zhu Y, Ulanova L, Begum N, Ding Y, Dau HA, et al. Matrix profile I: All
#'   pairs similarity joins for time series: A unifying view that includes motifs, discords and
#'   shapelets. Proc - IEEE Int Conf Data Mining, ICDM. 2017;1317–22.
#' @references * Zhu Y, Imamura M, Nikovski D, Keogh E. Matrix Profile VII: Time Series Chains: A
#'   New Primitive for Time Series Data Mining. Knowl Inf Syst. 2018 Jun 2;1–27.
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#' mp <- stamp.par(toy_data$data[1:200,1], window.size = 30, verbose = 0)
#' \dontrun{
#' ref.data <- toy_data$data[,1]
#' query.data <- toy_data$data[,2]
#' # self similarity
#' mp <- stamp.par(ref.data, window.size = 30, s.size = round(nrows(ref.data) * 0.1))
#' # join similarity
#' mp <- stamp.par(ref.data, query.data, window.size = 30, s.size = round(nrows(query.data) * 0.1))
#' }
#'
#' @import doSNOW foreach parallel
stamp.par <- function(..., window.size, exclusion.zone = 1 / 2, s.size = Inf, n.workers = 2, verbose = 2) {
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
    stop("Unknown type of data. Must be: a column matrix or a vector")
  }

  if (is.vector(query)) {
    query <- as.matrix(query)
  } else if (is.matrix(query)) {
    if (ncol(query) > nrow(query)) {
      query <- t(query)
    }
  } else {
    stop("Unknown type of query. Must be: a column matrix or a vector")
  }

  exclusion.zone <- floor(window.size * exclusion.zone)
  data.size <- nrow(data)
  query.size <- nrow(query)
  matrix.profile.size <- data.size - window.size + 1
  num.queries <- query.size - window.size + 1

  if (window.size > query.size / 2) {
    stop("Error: Time series is too short relative to desired subsequence length")
  }
  if (window.size < 4) {
    stop("Error: Subsequence length must be at least 4")
  }

  matrix.profile <- matrix(Inf, matrix.profile.size, 1)
  left.matrix.profile <- right.matrix.profile <- matrix.profile
  profile.index <- matrix(-1, matrix.profile.size, 1)
  left.profile.index <- right.profile.index <- profile.index

  ssize <- min(s.size, num.queries)
  order <- sample(1:num.queries, size = ssize)

  tictac <- Sys.time()

  cores <- min(max(2, n.workers), parallel::detectCores())

  cols <- min(num.queries, 100)

  lines <- 0:(ceiling(ssize / cols) - 1)
  if (verbose > 0) {
    pb <- utils::txtProgressBar(min = 0, max = max(lines), style = 3, width = 80)
  }
  cl <- parallel::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  on.exit(parallel::stopCluster(cl))
  if (verbose > 0) {
    on.exit(close(pb), TRUE)
  }
  if (verbose > 1) {
    on.exit(beep(sounds[[1]]), TRUE)
  }
  # anytime must return the result always
  on.exit(return(list(
    rmp = right.matrix.profile, rpi = right.profile.index,
    lmp = left.matrix.profile, lpi = left.profile.index,
    mp = matrix.profile, pi = profile.index
  )), TRUE)

  pre <- mass.pre(data, data.size, query, query.size, window.size = window.size)

  j <- NULL # CRAN NOTE fix
  `%dopar%` <- foreach::`%dopar%` # CRAN NOTE fix

  for (k in lines) {
    batch <- foreach(
      j = 1:cols,
      # .verbose = FALSE,
      .inorder = FALSE,
      .multicombine = TRUE,
      # .options.snow = opts,
      # .combine = combiner,
      # .errorhandling = 'remove',
      .export = "mass"
    ) %dopar% {
      res <- NULL

      index <- k * cols + j
      if (index <= ssize) {
        i <- order[index]
        nn <- mass(pre$data.fft, query[i:(i + window.size - 1)], data.size, window.size, pre$data.mean, pre$data.sd, pre$query.mean[i], pre$query.sd[i])
        distance.profile <- Re(sqrt(nn$distance.profile))

        if (exclusion.zone > 0) {
          distance.profile[max((i - exclusion.zone), 1):min((i + exclusion.zone), matrix.profile.size)] <- Inf
        }

        res <- list(dp = distance.profile, i = i)
      }

      res
    }

    for (i in 1:length(batch)) {
      curr <- batch[[i]]$i

      if (!is.null(curr)) {
        # left matrix.profile
        ind <- (batch[[i]]$dp[curr:matrix.profile.size] < left.matrix.profile[curr:matrix.profile.size])
        ind <- c(rep(FALSE, (curr - 1)), ind) # pad left
        left.matrix.profile[ind] <- batch[[i]]$dp[ind]
        left.profile.index[which(ind)] <- curr

        # right matrix.profile
        ind <- (batch[[i]]$dp[1:curr] < right.matrix.profile[1:curr])
        ind <- c(ind, rep(FALSE, matrix.profile.size - curr)) # pad right
        right.matrix.profile[ind] <- batch[[i]]$dp[ind]
        right.profile.index[which(ind)] <- curr

        # normal matrix.profile
        ind <- (batch[[i]]$dp < matrix.profile)
        matrix.profile[ind] <- batch[[i]]$dp[ind]
        profile.index[which(ind)] <- curr
      }
    }

    if (verbose > 0) {
      utils::setTxtProgressBar(pb, k)
    }
  }

  tictac <- Sys.time() - tictac

  if (verbose > 0) {
    message(sprintf("\nFinished in %.2f %s", tictac, units(tictac)))
  }

  # return() is at on.exit() function
}
