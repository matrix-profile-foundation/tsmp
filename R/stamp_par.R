#' Anytime univariate STAMP algorithm Parallel version
#'
#' No description
#'
#' No details
#'
#' @param ... a matrix or a vector, if a second time series is supplied it will be a join matrix profile
#' @param window.size size of the sliding window
#' @param exclusion.zone size of the exclusion zone, based on query size (default 1/2)
#' @param s.size for anytime algorithm, represents the size (in observations) the random calculation will occour (default Inf)
#'
#' @return The matrix profile and profile index
#' @export
#'
#' @family Stamp
#' @seealso [mstomp()]
#' @references Yeh CCM, Zhu Y, Ulanova L, Begum N, Ding Y, Dau HA, et al. Matrix profile I: All pairs similarity joins for time series: A unifying view that includes motifs, discords and shapelets. Proc - IEEE Int Conf Data Mining, ICDM. 2017;1317–22.
#' @references <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#' \dontrun{
#' mp <- stamp.par(data, window.size = 30)
#' mp <- stamp.par(ref.data, query.data, window.size = 30, s.size = round(nrows(ref.data) * 0.1))
#' }
#'
#' @import beepr doSNOW foreach parallel
stamp.par <- function(..., window.size, exclusion.zone = 1 / 2, s.size = Inf) {
  args <- list(...)
  data <- args[[1]]
  if (length(args) > 1)
    query <- args[[2]]
  else
    query <- data

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

  cores <- parallel::detectCores()
  cols <- min(data.size, 100)

  lines <- 0:(ceiling(ssize / cols) - 1)
  pb <- utils::txtProgressBar(min = 0, max = max(lines), style = 3)
  cl <- parallel::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  on.exit(parallel::stopCluster(cl))
  on.exit(close(pb), TRUE)
  on.exit(beepr::beep(), TRUE)
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
    batch <- foreach::foreach(
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
        distance.profile[max((i - exclusion.zone), 1):min((i + exclusion.zone - 1), matrix.profile.size)] <- Inf

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

    utils::setTxtProgressBar(pb, k)
  }

  # return() is at on.exit() function
}
