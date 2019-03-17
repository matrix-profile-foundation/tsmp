#' Real-time STAMP algorithm
#'
#' @param .mp a TSMP object of class `MatrixProfile`.
#' @param new_data new data to append to original data.
#'
#' @return Returns the input `.mp` updated with the new information.
#' @export
#'
#' @examples
stampi_update <- function(.mp, new_data) {
  new_data_size <- length(new_data)
  data_upd <- c(as.vector(.mp$data[[1]]), new_data)

  q1_idx <- (length(data_upd) - .mp$w + 1 - new_data_size + 1)

  mp_new <- c(.mp$mp, rep(Inf, new_data_size))
  pi_new <- c(.mp$pi, rep(-1, new_data_size))

  exclusion_zone <- (.mp$ez * .mp$w)

  for (i in seq_along(new_data)) {
    start_idx <- (q1_idx + i - 1)
    end_idx <- start_idx + .mp$w - 1
    query <- data_upd[start_idx:end_idx]

    nn <- dist_profile(data_upd, query)
    distance_profile <- abs(sqrt(nn$distance_profile))

    exc_st <- max(1, start_idx - exclusion_zone)
    exc_ed <- min(length(distance_profile), start_idx + exclusion_zone)
    distance_profile[exc_st:exc_ed] <- Inf

    upd_idxs <- distance_profile < mp_new

    pi_new[upd_idxs] <- start_idx
    mp_new[upd_idxs] <- distance_profile[upd_idxs]
    pi_new[start_idx] <- which.min(distance_profile)
    mp_new[start_idx] <- distance_profile[pi_new[start_idx]]
  }

  .mp$mp <- as.matrix(mp_new)
  .mp$pi <- as.matrix(pi_new)
  .mp$data[[1]] <- as.matrix(data_upd)

  return(.mp)
}

#' Real-time STOMP algorithm
#'
#' @param .mp a TSMP object of class `MatrixProfile`.
#' @param new_data new data to append to original data.
#'
#' @return Returns the input `.mp` updated with the new information.
#' @export
#'
#' @examples
stompi_update <- function(.mp, new_data) {
  new_data_size <- length(new_data)
  data_upd <- c(as.vector(.mp$data[[1]]), new_data)
  data_upd_size <- length(data_upd)

  q1_idx <- (data_upd_size - .mp$w + 1 - new_data_size + 1)

  mp_new <- c(.mp$mp, rep(Inf, new_data_size))
  pi_new <- c(.mp$pi, rep(-1, new_data_size))

  exclusion_zone <- (.mp$ez * .mp$w)

  if (new_data_size > 1) {
    rnn <- dist_profile(data_upd, data_upd[1:.mp$w])
    first_product <- rnn$last_product
    query_stats <- fast_avg_sd(data_upd[q1_idx:data_upd_size], .mp$w)
    drop_value <- 0
  }

  for (i in seq_len(new_data_size)) {
    start_idx <- (q1_idx + i - 1)
    end_idx <- start_idx + .mp$w - 1
    query <- data_upd[start_idx:end_idx]

    if (i == 1) {
      nn <- dist_profile(data_upd, query)
      distance_profile <- nn$distance_profile
      last_product <- nn$last_product
    } else {
      last_product[2:(data_upd_size - .mp$w + 1)] <- last_product[1:(data_upd_size - .mp$w)] -
        data_upd[1:(data_upd_size - .mp$w)] * drop_value +
        data_upd[(.mp$w + 1):data_upd_size] * query[.mp$w]
      last_product[1] <- first_product[start_idx]
      distance_profile <- 2 * (.mp$w - (last_product - .mp$w * nn$par$data_mean * query_stats$avg[i]) /
        (nn$par$data_sd * query_stats$sd[i]))
    }

    distance_profile <- abs(sqrt(distance_profile))
    drop_value <- query[1]

    exc_st <- max(1, start_idx - exclusion_zone)
    exc_ed <- min(length(distance_profile), start_idx + exclusion_zone)
    distance_profile[exc_st:exc_ed] <- Inf

    upd_idxs <- distance_profile < mp_new

    pi_new[upd_idxs] <- start_idx
    mp_new[upd_idxs] <- distance_profile[upd_idxs]
    pi_new[start_idx] <- which.min(distance_profile)
    mp_new[start_idx] <- distance_profile[pi_new[start_idx]]
  }

  .mp$mp <- as.matrix(mp_new)
  .mp$pi <- as.matrix(pi_new)
  .mp$data[[1]] <- as.matrix(data_upd)

  return(.mp)
}
