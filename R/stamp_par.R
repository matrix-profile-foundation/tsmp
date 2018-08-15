#' Anytime univariate STAMP algorithm Parallel version
#'
#' Computes the best so far Matrix Profile and Profile Index for Univariate Time Series.
#'
#' The Matrix Profile, has the potential to revolutionize time series data mining because of its generality, versatility, simplicity and scalability. In particular it has implications for time series motif discovery, time series joins, shapelet discovery (classification), density estimation, semantic segmentation, visualization, rule discovery, clustering etc.
#' The anytime STAMP computes the Matrix Profile and Profile Index in such manner that it can be stopped before its complete calculation and return the best so far results allowing ultra-fast approximate solutions.
#' `verbose` changes how much information is printed by this function; `0` means nothing, `1` means text, `2` means text and sound.
#'
#' @param ... a `matrix` or a `vector`. If a second time series is supplied it will be a join matrix profile.
#' @param window.size an `int`. Size of the sliding window.
#' @param exclusion.zone a `numeric`. Size of the exclusion zone, based on query size (default is `1/2`).
#' @param s.size a `numeric`. for anytime algorithm, represents the size (in observations) the random calculation will occour (default is `Inf`).
#' @param n.workers an `int`. Number of workers for parallel. (Default is `2`).
#' @param verbose an `int`. See details. (Default is `2`).
#'
#' @return Returns the matrix profile `mp` and profile index `pi`.
#' It also returns the left and right matrix profile `lmp`, `rmp` and profile index `lpi`, `rpi` that may be used to detect Time Series Chains (Yan Zhu 2018).
#' @export
#'
#' @family Stamp
#' @seealso [mstomp()], [mstomp.par()]
#' @references 1. Yeh CCM, Zhu Y, Ulanova L, Begum N, Ding Y, Dau HA, et al. Matrix profile I: All pairs similarity joins for time series: A unifying view that includes motifs, discords and shapelets. Proc - IEEE Int Conf Data Mining, ICDM. 2017;1317–22.
#' @references 2. Zhu Y, Imamura M, Nikovski D, Keogh E. Matrix Profile VII: Time Series Chains: A New Primitive for Time Series Data Mining. Knowl Inf Syst. 2018 Jun 2;1–27.
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#' Sys.sleep(1) # sometimes sleep is needed if you run parallel multiple times in a row
#' mp <- stamp.par(toy_data$data[1:200,1], window.size = 30, verbose = 0)
#' \dontrun{
#' mp <- stamp.par(ref.data, query.data, window.size = 30, s.size = round(nrows(ref.data) * 0.1))
#' }
#'
#' @import beepr doSNOW foreach parallel
stamp.par <- function(..., window.size, exclusion.zone = 1 / 2, s.size = Inf, n.workers = 2, verbose = 2) {
  args <- list(...)
  data <- args[[1]]
  if (length(args) > 1) {
    query <- args[[2]]
  } else {
    query <- data
  }

  if (is.vector(data)) {
    data <- as.matrix(data)
  }
  if (is.vector(query)) {
    query <- as.matrix(query)
  }

  if (is.matrix(data)) {
    if (ncol(data) > nrow(data)) {
      data <- t(data)
    }
  }

  if (is.matrix(query)) {
    if (ncol(query) > nrow(query)) {
      query <- t(query)
    }
  }

  exclusion.zone <- floor(window.size * exclusion.zone)
  data.size <- nrow(data)
  query.size <- nrow(data)
  matrix.profile.size <- data.size - window.size + 1

  matrix.profile <- matrix(Inf, matrix.profile.size, 1)
  left.matrix.profile <- right.matrix.profile <- matrix.profile
  profile.index <- matrix(-1, matrix.profile.size, 1)
  left.profile.index <- right.profile.index <- profile.index

  ssize <- min(s.size, matrix.profile.size)
  order <- sample(1:matrix.profile.size, size = ssize)

  cores <- min(max(2, n.workers), parallel::detectCores())

  cols <- min(data.size, 100)

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
    on.exit(beepr::beep(), TRUE)
  }
  # anytime must return the result always
  on.exit(return(list(
    rmp = right.matrix.profile, rpi = right.profile.index,
    lmp = left.matrix.profile, lpi = left.profile.index,
    mp = matrix.profile, pi = profile.index
  )), TRUE)

  pre <- mass.pre(data, data.size, query, query.size, window.size = window.size)

  tictac <- Sys.time()

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
      .export = c("mass")
    ) %dopar% {
      res <- NULL

      index <- k * cols + j
      if (index <= ssize) {
        i <- order[index]
        distance.profile <- Re(sqrt(mass(pre$data.fft, query[i:(i + window.size - 1)], data.size, window.size, pre$data.mean, pre$data.sd, pre$query.mean[i], pre$query.sd[i])$distance.profile))
        distance.profile[max((i - exclusion.zone), 1):min((i + exclusion.zone), matrix.profile.size)] <- Inf

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
