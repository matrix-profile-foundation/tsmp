#' Anytime univariate STAMP algorithm
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
#' mp <- stamp(data, window.size = 30)
#' mp <- stamp(ref.data, query.data, window.size = 30, s.size = round(nrows(ref.data) * 0.1))
#' }

stamp <- function(..., window.size, exclusion.zone = 1 / 2, s.size = Inf) {
  args <- list(...)
  data <- args[[1]]
  if (length(args) > 1)
    query <- args[[2]]
  else
    query <- data

  ## transform data list into matrix
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
  query.size <- nrow(data)
  matrix.profile.size <- data.size - window.size + 1

  matrix.profile <- matrix(Inf, matrix.profile.size, 1)
  left.matrix.profile <- right.matrix.profile <- matrix.profile
  profile.index <- matrix(-1, matrix.profile.size, 1)
  left.profile.index <- right.profile.index <- profile.index

  j <- 1
  ssize <- min(s.size, matrix.profile.size)
  order <- sample(1:matrix.profile.size, size = ssize)

  pb <- utils::txtProgressBar(min = 1, max = ssize, style = 3)

  on.exit(close(pb))
  on.exit(beepr::beep(), TRUE)
  # anytime must return the result always
  on.exit(return(list(
    rmp = as.matrix(right.matrix.profile), rpi = as.matrix(right.profile.index),
    lmp = as.matrix(left.matrix.profile), lpi = as.matrix(left.profile.index),
    mp = as.matrix(matrix.profile), pi = as.matrix(profile.index)
  )), TRUE)

  pre <- mass.pre(data, data.size, query, query.size, window.size)

  for (i in order) {
    j <- j + 1
    distance.profile <- Re(sqrt(mass(pre$data.fft, query[i:(i + window.size - 1), ], data.size, window.size, pre$data.mean, pre$data.sd, pre$query.mean[i], pre$query.sd[i])$distance.profile))
    distance.profile[max((i - exclusion.zone), 1):min((i + exclusion.zone - 1), matrix.profile.size)] <- Inf

    # anytime version
    # left matrix.profile
    ind <- (distance.profile[i:matrix.profile.size] < left.matrix.profile[i:matrix.profile.size])
    ind <- c(rep(FALSE, (i - 1)), ind) # pad left
    left.matrix.profile[ind] <- distance.profile[ind]
    left.profile.index[which(ind)] <- i

    # right matrix.profile
    ind <- (distance.profile[1:i] < right.matrix.profile[1:i])
    ind <- c(ind, rep(FALSE, matrix.profile.size - i)) # pad right
    right.matrix.profile[ind] <- distance.profile[ind]
    right.profile.index[which(ind)] <- i

    # normal matrix.profile
    ind <- (distance.profile < matrix.profile)
    matrix.profile[ind] <- distance.profile[ind]
    profile.index[which(ind)] <- i

    utils::setTxtProgressBar(pb, j)
  }

  # return() is at on.exit() function
}
