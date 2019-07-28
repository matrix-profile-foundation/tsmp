#' Real-time STOMP algorithm
#'
#' @param .mp a TSMP object of class `MatrixProfile`.
#' @param new_data new data to append to original data.
#' @param history_size an `int` or `FALSE`. (Default is `FALSE`). Keep only this amount of data in
#' the object. The value is for the data, not the matrix profile. Notice that the `lmp`and `lpi` will
#' be inconsistent when repeatedly updating limiting the history size and thus will affect
#' the `mp` and `pi`.
#'
#' @return Returns the input `.mp` updated with the new information.
#' @export
#'
#' @examples
#' mp <- tsmp(mp_toy_data$data[1:200, 1], window_size = 30, verbose = 0)
#' mpi <- stompi_update(mp, mp_toy_data$data[201:300, 1])
#' mp <- tsmp(mp_toy_data$data[1:300, 1], window_size = 30, verbose = 0)
#' all.equal(mp, mpi, check.attributes = FALSE)
stompi_update <- function(.mp, new_data, history_size = FALSE) {
  if (!is.null(attr(.mp, "join")) && attr(.mp, "join")) {
    stop("Update not implemented in Join similarity")
  }

  new_data_size <- length(new_data)

  if (new_data_size == 0) {
    stop("No new data")
  }

  data_upd <- c(as.vector(.mp$data[[1]]), new_data)
  data_upd_size <- length(data_upd)

  mp_new <- c(.mp$mp, rep(Inf, new_data_size))
  pi_new <- c(.mp$pi, rep(-Inf, new_data_size))
  lmp_new <- c(.mp$lmp, rep(Inf, new_data_size))
  lpi_new <- c(.mp$lpi, rep(-Inf, new_data_size))
  rmp_new <- c(.mp$rmp, rep(Inf, new_data_size))
  rpi_new <- c(.mp$rpi, rep(-Inf, new_data_size))

  q1_idx <- (data_upd_size - .mp$w + 1 - new_data_size + 1)

  mp_new_size <- length(mp_new)

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

    distance_profile[distance_profile < 0] <- 0
    distance_profile <- sqrt(distance_profile)
    drop_value <- query[1]

    exc_st <- max(1, start_idx - exclusion_zone)
    distance_profile[exc_st:mp_new_size] <- Inf

    upd_idxs <- distance_profile < mp_new
    pi_new[upd_idxs] <- start_idx
    mp_new[upd_idxs] <- distance_profile[upd_idxs]
    pi_new[start_idx] <- which.min(distance_profile)
    mp_new[start_idx] <- distance_profile[pi_new[start_idx]]

    # left matrix_profile
    upd_idxs <- (distance_profile[start_idx:mp_new_size] < lmp_new[start_idx:mp_new_size])
    upd_idxs <- c(rep(FALSE, (start_idx - 1)), upd_idxs) # pad left
    lpi_new[upd_idxs] <- start_idx
    lmp_new[upd_idxs] <- distance_profile[upd_idxs]
    lpi_new[start_idx] <- which.min(distance_profile)
    lmp_new[start_idx] <- distance_profile[lpi_new[start_idx]]

    # right matrix_profile
    upd_idxs <- (distance_profile[1:start_idx] < rmp_new[1:start_idx])
    upd_idxs <- c(upd_idxs, rep(FALSE, mp_new_size - start_idx)) # pad right
    rmp_new[upd_idxs] <- distance_profile[upd_idxs]
    rpi_new[upd_idxs] <- start_idx
  }

  if (history_size && (data_upd_size > history_size)) {
    data_upd <- utils::tail(data_upd, history_size)
    mp_new_size <- history_size - .mp$w + 1
    offset <- data_upd_size - history_size

    mp_new <- utils::tail(mp_new, mp_new_size)
    pi_new <- utils::tail(pi_new - offset, mp_new_size)
    lmp_new <- utils::tail(lmp_new, mp_new_size)
    lpi_new <- utils::tail(lpi_new - offset, mp_new_size)
    rmp_new <- utils::tail(rmp_new, mp_new_size)
    rpi_new <- utils::tail(rpi_new - offset, mp_new_size)

    if (is.null(attr(.mp, "offset"))) {
      attr(.mp, "offset") <- offset
    } else {
      attr(.mp, "offset") <- attr(.mp, "offset") + offset
    }

    # pi_new <- utils::tail(pi_new - attr(.mp, "offset"), mp_new_size)
  }

  .mp$mp <- as.matrix(mp_new)
  .mp$pi <- as.matrix(pi_new)
  .mp$lmp <- as.matrix(lmp_new)
  .mp$lpi <- as.matrix(lpi_new)
  .mp$rmp <- as.matrix(rmp_new)
  .mp$rpi <- as.matrix(rpi_new)
  .mp$data[[1]] <- as.matrix(data_upd)
  attr(.mp, "new_data") <- new_data_size

  # TODO: with tail or not (tail will recompute some things)
  #  if (history_size && (data_upd_size > history_size)) {
  #    return(utils::tail(.mp, history_size))
  #  } else {
  return(.mp)
  #  }
}
