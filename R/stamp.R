#' Anytime univariate STAMP algorithm
#'
#' No description
#'
#' No details
#'
#' @param data a matrix, with time series as a column
#' @param query.size size of the sliding window
#' @param exclusion.zone size of the exclusion zone, based on query size (default 1/2)
#' @param s.size for anytime algorithm, represents the size (in observations) the random calculation will occour (default Inf)
#'
#' @return The matrix profile and profile index
#' @export
#'
#' @seealso [mstomp()]
#' @references Yeh CCM, Zhu Y, Ulanova L, Begum N, Ding Y, Dau HA, et al. Matrix profile I: All pairs similarity joins for time series: A unifying view that includes motifs, discords and shapelets. Proc - IEEE Int Conf Data Mining, ICDM. 2017;1317–22.
#'
#' @examples
#' \dontrun{
#' mp <- stamp(data, 30)
#' mp <- stamp(data, 30, s.size = round(nrows(data) * 0.1)) # 10 percent of data
#' }

stamp <- function(data, query.size, exclusion.zone = 1 / 2, s.size = Inf) {
  exclusion.zone <- floor(query.size * exclusion.zone)
  data.size <- length(data)
  matrix.profile.size <- data.size - query.size + 1

  matrix.profile <- rep(Inf, matrix.profile.size)
  left.matrix.profile <- right.matrix.profile <- matrix.profile
  matrix.profile.index <- rep(-1, matrix.profile.size)
  left.matrix.profile.index <- right.matrix.profile.index <- matrix.profile.index

  j <- 1
  ssize <- min(s.size, matrix.profile.size)
  order <- sample(1:matrix.profile.size, size = ssize)

  pb <- txtProgressBar(min = 1, max = ssize, style = 3)

  on.exit(close(pb))
  on.exit(beep(), TRUE)
  # anytime must return the result always
  on.exit(return(list(rmp = right.matrix.profile, rmpi = right.matrix.profile.index,
                      lmp = left.matrix.profile, lmpi = left.matrix.profile.index,
                      mp = matrix.profile, mpi = matrix.profile.index)), TRUE)

  data.fft <- data
  data.fft[(data.size + 1):(query.size + data.size)] <- 0
  data.fft <- fft(data.fft) # precompute fft of ts
  data.mean <- MpMovAvg(data, query.size)
  data.sd <- MpMovsd(data, query.size) # precompute moving SD

  for (i in order) {
    j <- j + 1
    distance.profile <- MpDp(data.fft, data[i:(i + query.size - 1)], data.size, query.size, data.mean, data.sd, data.mean[i], data.sd[i])
    distance.profile[max((i - exclusion.zone), 1):min((i + exclusion.zone - 1), matrix.profile.size)] <- Inf

    # anytime version
    # left matrix.profile
    ind <- (distance.profile[i:matrix.profile.size] < left.matrix.profile[i:matrix.profile.size])
    ind <- c(rep(FALSE, (i - 1)), ind) # pad left
    left.matrix.profile[ind] <- distance.profile[ind]
    left.matrix.profile.index[which(ind)] <- i

    # right matrix.profile
    ind <- (distance.profile[1:i] < right.matrix.profile[1:i])
    ind <- c(ind, rep(FALSE, matrix.profile.size - i)) # pad right
    right.matrix.profile[ind] <- distance.profile[ind]
    right.matrix.profile.index[which(ind)] <- i

    # normal matrix.profile
    ind <- (distance.profile < matrix.profile)
    matrix.profile[ind] <- distance.profile[ind]
    matrix.profile.index[which(ind)] <- i

    setTxtProgressBar(pb, j)
  }

  # return() is at on.exit() function
}
