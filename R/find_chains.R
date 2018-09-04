#' Find Time Series Chains
#'
#' Time Series Chains is a new primitive for time series data mining.
#'
#' @param matrices a result from STAMP or STOMP algorithms
#'
#' @return Returns `chains`, a `list` of chains founded with more than 2 patterns and `best` with
#'   the best one.
#' @export
#' @references * Zhu Y, Imamura M, Nikovski D, Keogh E. Matrix Profile VII: Time Series Chains: A
#'   New Primitive for Time Series Data Mining. Knowl Inf Syst. 2018 Jun 2;1–27.
#' @references Website: <https://sites.google.com/site/timeserieschain/>
#' @examples
#' w <- 50
#' data <- mp_gait_data
#' mp <- stamp(data, window.size = w, exclusion.zone = 1/4, verbose = 0)
#' find.chains(mp)
#'
find_chains <- function(matrices) {
  size <- length(matrices$rpi)
  chain_length <- rep(1, size)
  chain_set <- list()

  k <- 1

  for (i in 1:size) {
    if (chain_length[i] == 1) {
      j <- i
      chain <- j

      while (matrices$rpi[j] > 0 && matrices$lpi[matrices$rpi[j]] == j) {
        j <- matrices$rpi[j]
        chain_length[j] <- -1
        chain_length[i] <- chain_length[i] + 1
        chain <- c(chain, j)
      }

      if (length(chain) > 2) {
        chain_set[[k]] <- chain
        k <- k + 1
      }
    }
  }

  l <- max(chain_length)

  best_chain <- NULL
  mean <- Inf
  for (i in 1:length(chain_set)) {
    if (length(chain_set[[i]]) == l) {
      n <- mean(matrices$rmp[chain_set[[i]]])
      if (n < mean) {
        mean <- n
        best_chain <- chain_set[[i]]
      }
    }
  }

  return(list(chains = chain_set, best = best_chain))
}
