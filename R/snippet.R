## Anonymous-Author information blinded for review
## This is the source code for the ICDM paper "Time Series Snippets: A New Primitive
## for Time Series Data Mining". For more details, please refer to the
## supporting website: https://sites[['google']][['com']]/site/snippetfinder

## input:
## data : Time Series
## N : number of snippets the user wishes to find
## sub : the length of snippet
## per : the MPdist subsequence length percentage.
## For example, if per is 100 then the MPdist subsequence is the same as sub.

## output:
## snippet : a list of N snippets
## fraction : fraction of each snippet
## snippetidx : the location of each
## snippet within time sereis


#' Title
#'
#' @param data
#' @param N
#' @param sub
#' @param per
#'
#' @return
#' @export
#'
#' @examples
snippetfinder <- function(data, N, sub, per) {

  # currently is about 5.1x slower than MATLAB, but is working

  ## check input
  if (length(data) < 2 * sub) {
    stop("Error: Time series is too short relative to desired snippet Length")
  }

  ## initialization
  distances <- NULL
  indexes <- NULL
  snippet <- NULL
  minis <- Inf
  fraction <- NULL
  snippetidx <- NULL
  distancesSnipp <- NULL

  len <- sub * ceiling(length(data) / sub) - length(data)
  ## padding zeros to the end of the time series.
  data <- c(data, rep(len, 0))

  ## compute all the profiles
  for (i in seq.int(1, length(data) - sub, sub)) {
    indexes <- c(indexes, i)
    distance <- mpdist(data, data[i:(i + sub)], round(sub * per / 100))
    distances <- rbind(distances, distance)
  }

  ## calculating the Nth snippets
  for (n in 1:N) {
    minims <- Inf
    for (i in seq_len(nrow(distances))) {
      if (minims > sum(pmin(distances[i, ], minis))) {
        minims <- sum(pmin(distances[i, ], minis))
        index <- i
      }
    }
    minis <- pmin(distances[index, ], minis)
    snippet <- c(snippet, data[indexes[index]:(indexes[index] + sub)])
    snippetidx <- c(snippetidx, indexes[index])

    distance <- mpdist(data, data[indexes[index]:(indexes[index] + sub)], round(sub * per / 100))
    distancesSnipp <- rbind(distancesSnipp, distance)

    graphics::plot(distance, type = "l", ylab = "distance",
                   main = paste0('MPdist-', n, ' location-', indexes[index]))

    graphics::plot(data[indexes[index]:(indexes[index] + sub)], type = "l", ylab = "data",
                   main = paste0('snippet-', n, ' location-', indexes[index]))


    # box off;xlim([0 length(data(indexes(index):indexes(index)+sub))])
  }
  ## Calculating the fraction of each snippet
  totalmin <- do.call(pmin, as.data.frame(t(distancesSnipp)))

  for (i in 1:N) {
    a <- (distancesSnipp[i, ] <= totalmin)
    fraction <- c(fraction, sum(a) / (length(data) - sub))
    # # # #     just in case we have the same value for both
    totalmin[a] <- totalmin[a] - 1
  }

  ## Calculating the horizantal regime bar for 2 snippets
  if (N == 2) {
    totalmin <- do.call(pmin, as.data.frame(t(distancesSnipp)))

    a <- (distancesSnipp[1, ] <= totalmin)

    for (i in seq.int(1, length(a) - sub, sub)) {
      if (sum(a[i:(i + sub - 1)]) > (1 / 2 * sub)) {
        a[i:(i + sub - 1)] <- TRUE
      } else {
        a[i:(i + sub - 1)] <- FALSE
      }
    }

    # # # # # # # # # # # # # # # # #
    c <- matrix(0, nrow = length(a), ncol = 3)
    j <- 1
    i <- 1

    while (i <= length(a)) {
      b <- which(a[(i + 1):length(a)] != a[i])[1]
      if (is.na(b)) {
        break
      } else {
        c[j, ] <- c(i, i + b - 1, a[i])
        j <- j + 1
        i <- i + b
      }
    }
    if (i <= length(a)) {
      c[j, ] <- c(i, length(a), a[i])
    }
    c <- c[1:j, , drop = FALSE]

    graphics::plot(0.5, 0.5, ylab = "", xlab = "Index",
      type = "n", main = "Horizontal regime bar", xlim = c(min(c[, 1]), max(c[, 2])), ylim = c(0, 2)
    )

    for (i in 1:nrow(c)) {
      d <- (c[i, 1]:c[i, 2])
      if (c[i, 3] == 0) {
        graphics::lines(d, rep(1, length(d)), col = "blue")
      } else {
        graphics::lines(d, rep(1, length(d)), col = "red")
      }
    }
  }

  return(list(fraction = fraction, snippet = snippet, snippetidx = snippetidx))
}
