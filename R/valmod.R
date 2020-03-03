#' Variable Length Motif Discovery
#'
#' Computes the Matrix Profile and Profile Index for a range of query window sizes
#'
#' @details
#' This algorithm uses an exact algorithm based on a novel lower bounding technique, which is
#' specifically designed for the motif discovery problem. `verbose` changes how much information
#' is printed by this function; `0` means nothing, `1` means text, `2` adds the progress bar,
#' `3` adds the finish sound. `exclusion_zone` is used to avoid  trivial matches; if a query data
#' is provided (join similarity), this parameter is ignored.
#'
#' Paper that implements `skimp()` suggests that window_max / window_min > than 1.24 begins to
#' weakening pruning in `valmod()`.
#'
#' @param \dots a `matrix` or a `vector`. If a second time series is supplied it will be a join matrix
#'   profile.
#' @param window_min an `int`. Minimum size of the sliding window.
#' @param window_max an `int`. Maximum size of the sliding window.
#' @param heap_size an `int`. (Default is `50`). Size of the distance profile heap buffer
#' @param exclusion_zone a `numeric`. Size of the exclusion zone, based on window size (default is
#'   `1/2`). See details.
#' @param lb a `logical`. (Default is `TRUE`). If `FALSE` all window sizes will be calculated using
#' STOMP instead of pruning. This is just for academic purposes.
#' @param verbose an `int`. See details. (Default is `2`).
#'
#' @return Returns a `Valmod` object, a `list` with the matrix profile `mp`, profile index `pi`
#'   left and right matrix profile `lmp`, `rmp` and profile index `lpi`, `rpi`, best window size `w`
#'   for each index and exclusion zone `ez`. Additionally: `evolution_motif` the best motif distance
#'   per window size, and non-length normalized versions of `mp`, `pi` and `w`: `mpnn`, `pinn` and `wnn`.
#'
#' @export
#'
#' @family matrix profile computations
#'
#' @references * Linardi M, Zhu Y, Palpanas T, Keogh E. VALMOD: A Suite for Easy and Exact Detection
#'  of Variable Length Motifs in Data Series. In: Proceedings of the 2018 International Conference
#'   on Management of Data - SIGMOD '18. New York, New York, USA: ACM Press; 2018. p. 1757-60.
#'
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#' mp <- valmod(mp_toy_data$data[1:200, 1], window_min = 30, window_max = 40, verbose = 0)
#' \dontrun{
#' ref_data <- mp_toy_data$data[, 1]
#' query_data <- mp_toy_data$data[, 2]
#' # self similarity
#' mp <- valmod(ref_data, window_min = 30, window_max = 40)
#' # join similarity
#' mp <- valmod(ref_data, query_data, window_min = 30, window_max = 40)
#' }
#'
valmod <- function(..., window_min, window_max, heap_size = 50, exclusion_zone = getOption("tsmp.exclusion_zone", 1 / 2), lb = TRUE, verbose = getOption("tsmp.verbose", 2)) {
  argv <- list(...)
  argc <- length(argv)
  data <- argv[[1]]
  if (argc > 1 && !is.null(argv[[2]])) {
    query <- argv[[2]]
    exclusion_zone <- 0 # don't use exclusion zone for joins
    join <- TRUE
  } else {
    query <- data
    join <- FALSE
  }

  # transform data into matrix
  if (is.vector(data)) {
    data <- as.matrix(data)
  } else if (is.matrix(data)) {
    if (ncol(data) > nrow(data)) {
      data <- t(data)
    }
  } else {
    stop("Unknown type of data. Must be: a column matrix or a vector.", call. = FALSE)
  }

  if (is.vector(query)) {
    query <- as.matrix(query)
  } else if (is.matrix(query)) {
    if (ncol(query) > nrow(query)) {
      query <- t(query)
    }
  } else {
    stop("Unknown type of query. Must be: a column matrix or a vector.", call. = FALSE)
  }

  ez <- exclusion_zone # store original
  data_size <- nrow(data)
  query_size <- nrow(query)

  if (data_size != query_size) {
    stop("Join similarity for different sizes not implemented yet")
  }

  range_size <- window_max - window_min + 1
  max_profile_size <- data_size - window_min + 1

  if (window_min > query_size / 2) {
    stop("Time series is too short relative to desired window size.", call. = FALSE)
  }
  if (window_min < 4) {
    stop("`window_size` must be at least 4.", call. = FALSE)
  }

  # check skip position
  data_skip <- (is.na(data) | is.infinite(data))
  query_skip <- (is.na(query) | is.infinite(query))

  data[data_skip] <- 0
  query[query_skip] <- 0

  valmp <- list(
    matrix_profile = matrix(Inf, max_profile_size, 1),
    profile_index = matrix(-1, max_profile_size, 1),
    length_profile = matrix(-1, max_profile_size, 1),
    distances_evolution_motif = matrix(Inf, range_size, 1),
    matrix_profile_non_length_normalized = matrix(Inf, max_profile_size, 1),
    profile_index_non_length_normalized = matrix(-1, max_profile_size, 1),
    length_profile_non_length_normalized = matrix(-1, max_profile_size, 1)
  )

  list_motifs_profile <- array(0,
    dim = c(max_profile_size, 10, heap_size),
    dimnames = list(NULL, vars = c(
      "distances", "query_sd", "sum_query", "sum_data",
      "sqrsum_query", "sqrsum_data", "lb_distances",
      "index_query", "indexes_data", "dps"
    ), NULL)
  )

  motifs_per_size <- -1
  min_number_motifs_found <- -1
  max_number_motifs_found <- -1
  list_lb_min_non_valid <- vector(mode = "numeric", length = max_profile_size)
  index_lb_min_non_valid <- vector(mode = "numeric", length = max_profile_size)

  tictac <- Sys.time()

  # Pre computations
  data_stats <- vector(mode = "list", length = range_size)
  query_stats <- vector(mode = "list", length = range_size)
  for (i in seq_len(range_size)) {
    window_size <- window_min + i - 1
    data_stats[[i]] <- fast_avg_sd(data, window_size)
    query_stats[[i]] <- fast_avg_sd(query, window_size)
  }

  if (verbose > 1) {
    pbp <- progress::progress_bar$new(
      format = ":what [:bar] :percent at :tick_rate it/s, elapsed: :elapsed, eta: :eta",
      clear = FALSE, total = range_size, width = 80
    )
  }

  for (window_size in seq(window_min, window_max, 1)) {
    # ============== Valmod loop ===============

    exclusion_zone <- round(window_size * ez + vars()$eps)
    matrix_profile_size <- data_size - window_size + 1
    num_queries <- query_size - window_size + 1
    offset <- window_size - window_min

    skip_location <- rep(FALSE, matrix_profile_size)

    for (i in 1:matrix_profile_size) {
      if (any(data_skip[i:(i + window_size - 1)])) {
        skip_location[i] <- TRUE
      }
    }

    if (offset == 0 || motifs_per_size == 0 || lb == FALSE) {
      # ==== STOMP ====

      if (verbose == 1) {
        message("=== STOMP ===")
      }

      if (verbose > 1) {
        pbp$tick(0, tokens = list(what = "STOMP  "))
        pb <- progress::progress_bar$new(
          format = "STOMP [:bar] :percent at :tick_rate it/s, elapsed: :elapsed, eta: :eta",
          clear = FALSE, total = num_queries, width = 80
        )
      }

      if (verbose > 2) {
        on.exit(beep(sounds[[1]]), TRUE)
      }

      first_product <- matrix(0, num_queries, 1)

      # forward
      nn <- dist_profile(data, query, window_size = window_size)

      # reverse
      # This is needed to handle with the join similarity.
      rnn <- dist_profile(query, data, window_size = window_size)

      first_product[, 1] <- rnn$last_product

      matrix_profile <- matrix(Inf, matrix_profile_size, 1)
      profile_index <- matrix(-Inf, matrix_profile_size, 1)
      if (join) {
        # no RMP and LMP for joins
        left_matrix_profile <- right_matrix_profile <- NULL
        left_profile_index <- right_profile_index <- NULL
      } else {
        left_matrix_profile <- right_matrix_profile <- matrix_profile
        left_profile_index <- right_profile_index <- profile_index
      }
      distance_profile <- matrix(0, matrix_profile_size, 1)
      last_product <- matrix(0, matrix_profile_size, 1)
      drop_value <- matrix(0, 1, 1)

      for (i in seq_len(num_queries)) {
        # compute the distance profile
        query_window <- as.matrix(query[i:(i + window_size - 1), 1])

        if (i == 1) {
          distance_profile[, 1] <- nn$distance_profile
          last_product[, 1] <- nn$last_product
        } else {
          last_product[2:(data_size - window_size + 1), 1] <- last_product[1:(data_size - window_size), 1] -
            data[1:(data_size - window_size), 1] * drop_value +
            data[(window_size + 1):data_size, 1] * query_window[window_size, 1]

          last_product[1, 1] <- first_product[i, 1]
          distance_profile <- 2 * (window_size - (last_product - window_size * nn$par$data_mean * nn$par$query_mean[i]) /
            (nn$par$data_sd * nn$par$query_sd[i]))
        }

        distance_profile[distance_profile < 0] <- 0

        new_lb_profile <- ((last_product / window_size) - (nn$par$query_mean[i] * nn$par$data_mean)) / (nn$par$query_sd[i] * nn$par$data_sd)
        new_lb_profile[new_lb_profile < 0] <- 0
        new_lb_profile_len <- length(new_lb_profile)

        if (new_lb_profile_len != matrix_profile_size) {
          stop("new_lb_profile_len != matrix_profile_size")
        }

        lb_profile <- new_lb_profile
        lb_profile[is.na(lb_profile)] <- Inf
        lb_idx <- (lb_profile > 0)
        lb_profile[lb_idx] <- window_size * (1 - (lb_profile[lb_idx]^2))
        lb_profile[!lb_idx] <- window_size

        drop_value <- query_window[1, 1]

        # apply exclusion zone
        if (exclusion_zone > 0) {
          exc_st <- max(1, i - exclusion_zone)
          exc_ed <- min(matrix_profile_size, i + exclusion_zone)

          distance_profile[exc_st:exc_ed] <- Inf
          lb_profile[exc_st:exc_ed] <- Inf
        }

        distance_profile[nn$par$data_sd < vars()$eps] <- Inf
        lb_profile[nn$par$data_sd < vars()$eps] <- Inf
        if (skip_location[i] || any(nn$par$query_sd[i] < vars()$eps)) {
          distance_profile[] <- Inf
          lb_profile[] <- Inf
        }

        distance_profile[skip_location] <- Inf
        lb_profile[skip_location] <- Inf

        # ==== Store best 'p' LB DPs  ====

        lb_idxs <- sort.int(lb_profile, index.return = TRUE)$ix[1:heap_size]

        if (any(lb_idxs > matrix_profile_size)) {
          stop("lb_idxs > matrix_profile_size")
        }

        list_motifs_profile[i, "distances", ] <- distance_profile[lb_idxs]
        list_motifs_profile[i, "query_sd", ] <- nn$par$query_sd[i]
        list_motifs_profile[i, "sum_query", ] <- query_stats[[offset + 1]]$sum[i]
        list_motifs_profile[i, "sum_data", ] <- data_stats[[offset + 1]]$sum[lb_idxs]
        list_motifs_profile[i, "sqrsum_query", ] <- query_stats[[offset + 1]]$sqrsum[i]
        list_motifs_profile[i, "sqrsum_data", ] <- data_stats[[offset + 1]]$sqrsum[lb_idxs]
        list_motifs_profile[i, "lb_distances", ] <- lb_profile[lb_idxs]
        list_motifs_profile[i, "index_query", ] <- i
        list_motifs_profile[i, "indexes_data", ] <- lb_idxs
        list_motifs_profile[i, "dps", ] <- last_product[lb_idxs]

        distance_profile <- sqrt(distance_profile)

        ind <- (distance_profile < matrix_profile)
        matrix_profile[is.na(ind)] <- NA
        ind <- ind & !is.na(ind)
        matrix_profile[ind] <- distance_profile[ind]
        profile_index[which(ind)] <- i

        if (verbose > 1) {
          pb$tick()
        }
      }

      # ==== Update VALMP  ====
      ind <- (matrix_profile < valmp$matrix_profile_non_length_normalized[1:matrix_profile_size])
      valmp$matrix_profile_non_length_normalized[which(ind)] <- matrix_profile[which(ind)]
      valmp$profile_index_non_length_normalized[which(ind)] <- profile_index[which(ind)]
      valmp$length_profile_non_length_normalized[ind] <- window_size
      ## normalize distance
      normalized_distance <- matrix_profile * sqrt(1.0 / window_size)
      ind <- (normalized_distance < valmp$matrix_profile[1:matrix_profile_size])
      valmp$matrix_profile[which(ind)] <- normalized_distance[which(ind)]
      valmp$profile_index[which(ind)] <- profile_index[which(ind)]
      valmp$length_profile[ind] <- window_size

      valmp$distances_evolution_motif[offset + 1, ] <- min(matrix_profile)^2 * sqrt(1.0 / window_size)

      motifs_per_size <- -1

      if (verbose > 1) {
        pbp$tick()

        if (!pb$finished) {
          pb$terminate()
        }
      }

      if (verbose == 1) {
        message("=== PRUNING ===")
      }
    } else {
      #### LB Pruning ####

      if (verbose > 1) {
        pbp$tick(0, tokens = list(what = "Pruning "))
      }

      if (offset == 0) {
        warning("OFFSET ZERO?")
      }
      # offset never gets zero in this block, why check offset == 0 ?
      motifs_per_size <- 0

      list_valid_entries <- vector(mode = "numeric", length = matrix_profile_size)
      index_valid_entries <- vector(mode = "numeric", length = matrix_profile_size)
      non_valid_smaller <- -1
      min_abs_true_dist <- -1

      i_v <- seq_len(matrix_profile_size)
      # --- First Outer Loop ----

      new_index_query_v <- list_motifs_profile[i_v, "index_query", 1] + window_size - 1

      i_v <- i_v[new_index_query_v <= query_size]
      new_index_query_v <- new_index_query_v[new_index_query_v <= query_size]

      max_lb <- list_motifs_profile[i_v, "lb_distances", heap_size]
      max_query_sd <- list_motifs_profile[i_v, "query_sd", heap_size]
      curr_query_mean <- query_stats[[offset + 1]]$avg[i_v]
      curr_query_sd <- query_stats[[offset + 1]]$sd[i_v]
      lower_bound <- max_lb * max_query_sd^2 / curr_query_sd^2

      min_entry_idx <- NULL
      min_entry_lb <- -1

      list_motifs_profile[i_v, "sum_query", ] <- list_motifs_profile[i_v, "sum_query", ] + query[new_index_query_v]
      list_motifs_profile[i_v, "sqrsum_query", ] <- list_motifs_profile[i_v, "sqrsum_query", ] + (query[new_index_query_v] * query[new_index_query_v])
      query_mean_v <- matrix(rep(curr_query_mean, heap_size), ncol = heap_size)
      query_sd_v <- matrix(rep(curr_query_sd, heap_size), ncol = heap_size)

      # --- First Inner Loop ----
      j_v <- seq(heap_size, by = -1)

      # apply exclusion zone
      # !ezx_v contains all trival matches to recompute the DP
      ezx_v <- (((list_motifs_profile[i_v, "indexes_data", j_v]) < (list_motifs_profile[i_v, "index_query", j_v] - exclusion_zone)) |
        ((list_motifs_profile[i_v, "indexes_data", j_v]) > (list_motifs_profile[i_v, "index_query", j_v] + exclusion_zone)))

      ez_v <- ezx_v & (list_motifs_profile[i_v, "indexes_data", j_v] + window_size - 1 <= data_size)

      new_index_data_v <- list_motifs_profile[i_v, "indexes_data", j_v] + window_size - 1

      # --- Compute true distance entry heapentry ----
      # TODO: check NA from data
      list_motifs_profile[i_v, "dps", j_v][ez_v] <- (list_motifs_profile[i_v, "dps", j_v] + matrix(query[new_index_query_v] * data[new_index_data_v], ncol = ncol(new_index_data_v)))[ez_v]
      list_motifs_profile[i_v, "sum_data", j_v][ez_v] <- (list_motifs_profile[i_v, "sum_data", j_v] + data[new_index_data_v])[ez_v]
      list_motifs_profile[i_v, "sqrsum_data", j_v][ez_v] <- (list_motifs_profile[i_v, "sqrsum_data", j_v] + (data[new_index_data_v] * data[new_index_data_v]))[ez_v]

      data_mean_v <- list_motifs_profile[i_v, "sum_data", ] / window_size
      data_sd_v <- (list_motifs_profile[i_v, "sqrsum_data", ] / window_size) - (data_mean_v * data_mean_v)


      data_sd_v[data_sd_v < 0] <- 0
      data_sd_v <- sqrt(data_sd_v)

      dist_v <- ifelse(ez_v,
        ((2 * window_size) * (1 - ((list_motifs_profile[i_v, "dps", j_v] -
          (window_size * query_mean_v[i_v, j_v] * data_mean_v[i_v, j_v])) /
          (window_size * query_sd_v[i_v, j_v] * data_sd_v[i_v, j_v])))),
        list_motifs_profile[i_v, "distances", j_v]
      )

      dist_v[dist_v < 0] <- 0

      list_motifs_profile[i_v, "distances", j_v] <- dist_v

      # ---- Define min_entry ----
      distances <- list_motifs_profile[i_v, "distances", j_v]
      distances[!ez_v] <- Inf
      min_entry_idx <- j_v[apply(distances, 1, which.min)]
      lbdistances <- list_motifs_profile[i_v, "lb_distances", j_v]
      lbdistances[!ez_v] <- 0
      lbdistances <- lbdistances * list_motifs_profile[i_v, "query_sd", j_v]^2 / curr_query_sd^2
      min_entry_lb <- do.call(pmax, lapply(seq_len(ncol(lbdistances)), function(i) lbdistances[, i]))


      #### Pruning is effective global min found ####
      # update the min OOOOONLYYYY IF it is a correct min of the distance profile
      min_distances <- list_motifs_profile[i_v, "distances", , drop = FALSE][cbind(seq_along(min_entry_idx), 1, min_entry_idx)]
      valid_entries <- min_distances < lower_bound

      if (any(valid_entries)) {
        min_abs_true_dist <- min(min_distances)
        list_valid_entries <- i_v[valid_entries]
        index_valid_entries <- min_entry_idx[valid_entries]
      }

      #### The distance profile has elements but nothing to say, store the min max LB ####
      if (any(!valid_entries)) {
        # store the non valid entry with the smallest LB.
        # (I am sure that this guy correctly lowerbounds the MP  )
        non_valid_entries <- (!valid_entries) & (min_entry_lb >= 0)
        list_lb_min_non_valid <- min_entry_lb[non_valid_entries]
        index_lb_min_non_valid <- i_v[non_valid_entries]
        non_valid_smaller <- min(list_lb_min_non_valid)
      }

      all_trivial <- which(apply(ezx_v, 1, sum) == 0)

      if (length(all_trivial) > 0) {
        if (length(i_v[!(i_v[all_trivial] %in% i_v[!valid_entries])]) > 0) {
          warning("Warning: Some trivial matches might not have been recomputed")
        }
      }

      # check to see how many valid motifs (smallest value of the Matrix Profile) I have among
      # the valid entries.

      best_motif <- NULL

      # --- Second Outer Loop ----
      m <- (list_motifs_profile[list_valid_entries, "distances", , drop = FALSE][cbind(seq_along(index_valid_entries), 1, index_valid_entries)] < non_valid_smaller)
      motifs_per_size <- sum(m)

      if (motifs_per_size > 0) {
        i <- list_valid_entries[m]
        j <- index_valid_entries[m]
        distances <- list_motifs_profile[i, "distances", , drop = FALSE][cbind(seq_along(j), 1, j)]

        # UPDATE VALMAP_t!
        real_distance <- sqrt(distances)
        normalized_distance <- sqrt(distances) * sqrt(1.0 / window_size)

        n <- (normalized_distance < valmp$matrix_profile[i])

        if (any(n)) {
          entries <- i[n]
          indexes <- list_motifs_profile[entries, "indexes_data", , drop = FALSE][cbind(seq_along(j[n]), 1, j[n])]

          valmp$matrix_profile[entries] <- normalized_distance[n]
          valmp$profile_index[entries] <- indexes
          valmp$length_profile[entries] <- window_size
        }

        r <- (real_distance < valmp$matrix_profile_non_length_normalized[i])

        if (any(r)) {
          entries <- i[r]
          indexes <- list_motifs_profile[entries, "indexes_data", , drop = FALSE][cbind(seq_along(j[r]), 1, j[r])]

          valmp$matrix_profile_non_length_normalized[entries] <- real_distance[r]
          valmp$profile_index_non_length_normalized[entries] <- indexes
          valmp$length_profile_non_length_normalized[entries] <- window_size
        }

        # TODO: check if offset MUST be reduced when re-STOMP

        best_motif <- min(distances)
      }

      if (is.null(best_motif)) {
        # matrix profile may not be containing the minimum
        # Here I do not want to call STOMP just update the distance profiles which have the max
        # LB smaller than the minimum motif pair

        non_valid_idxs <- (list_lb_min_non_valid <= min_abs_true_dist)

        if (any(non_valid_idxs)) {
          indexes_to_update <- sort(index_lb_min_non_valid[non_valid_idxs])

          if (verbose > 1) {
            pbp$tick(0, tokens = list(what = sprintf("DPs  %s  ", length(indexes_to_update))))
          }

          # check for sequences where we can use STOMP instead of MASS
          idif <- c(1, diff(indexes_to_update))
          idxs <- c(1, which(idif > 1), length(indexes_to_update) + 1)

          sequences <- list()

          for (i in seq_len(length(idxs) - 1)) {
            seq <- indexes_to_update[idxs[i]:(idxs[i + 1] - 1)]
            sequences[[i]] <- seq
          }

          nn <- dist_profile(data, query, window_size = window_size)
          rnn <- dist_profile(query, data, window_size = window_size)

          first_product <- rnn$last_product

          for (j in seq_along(sequences)) {
            seq <- sequences[[j]]

            for (i in seq_along(seq)) {
              query_window <- query[seq[i]:(seq[i] + window_size - 1)]

              if (i == 1) {
                if (verbose > 1) {
                  pbp$tick(0, tokens = list(what = sprintf("MASS    ")))
                }

                nni <- dist_profile(data, query, nn, index = seq[i])
                distance_profile <- nni$distance_profile
                last_product <- nni$last_product
              } else {
                if (verbose > 1) {
                  pbp$tick(0, tokens = list(what = sprintf("STOMP %s ", i)))
                }

                last_product[2:(data_size - window_size + 1)] <- last_product[1:(data_size - window_size)] -
                  data[1:(data_size - window_size)] * drop_value +
                  data[(window_size + 1):data_size] * query_window[window_size]

                last_product[1] <- first_product[seq[i]]
                distance_profile <- 2 * (window_size - (last_product - window_size * nn$par$data_mean * nn$par$query_mean[seq[i]]) /
                  (nn$par$data_sd * nn$par$query_sd[seq[i]]))
              }

              drop_value <- query_window[1]

              distance_profile[distance_profile < 0] <- 0

              new_lb_profile <- ((last_product / window_size) - (nn$par$query_mean[seq[i]] *
                nn$par$data_mean)) / (nn$par$query_sd[seq[i]] * nn$par$data_sd)
              new_lb_profile[new_lb_profile < 0] <- 0
              new_lb_profile_len <- length(new_lb_profile)

              if (new_lb_profile_len != matrix_profile_size) {
                stop("new_lb_profile_len != matrix_profile_size")
              }

              lb_profile <- new_lb_profile
              lb_idx <- (lb_profile > 0)
              lb_profile[lb_idx] <- window_size * (1 - (lb_profile[lb_idx]^2))
              lb_profile[!lb_idx] <- window_size

              # apply exclusion zone
              if (exclusion_zone > 0) {
                exc_st <- max(1, seq[i] - exclusion_zone)
                exc_ed <- min(matrix_profile_size, seq[i] + exclusion_zone)

                distance_profile[exc_st:exc_ed] <- Inf
                lb_profile[exc_st:exc_ed] <- Inf
              }

              distance_profile[nn$par$data_sd < vars()$eps] <- Inf
              lb_profile[nn$par$data_sd < vars()$eps] <- Inf
              if (skip_location[seq[i]] || any(nn$par$query_sd[seq[i]] < vars()$eps)) {
                distance_profile[] <- Inf
                lb_profile[] <- Inf
              }

              distance_profile[skip_location] <- Inf
              lb_profile[skip_location] <- Inf

              lb_idxs <- sort.int(lb_profile, index.return = TRUE)$ix[1:heap_size]

              if (any(lb_idxs > matrix_profile_size)) {
                stop("lb_idxs > matrix_profile_size")
              }

              list_motifs_profile[seq[i], "distances", ] <- distance_profile[lb_idxs]
              list_motifs_profile[seq[i], "query_sd", ] <- query_sd_v[seq[i], ]
              list_motifs_profile[seq[i], "sum_query", ] <- query_stats[[offset + 1]]$sum[seq[i]]
              list_motifs_profile[seq[i], "sum_data", ] <- data_stats[[offset + 1]]$sum[lb_idxs]
              list_motifs_profile[seq[i], "sqrsum_query", ] <- query_stats[[offset + 1]]$sqrsum[seq[i]]
              list_motifs_profile[seq[i], "sqrsum_data", ] <- data_stats[[offset + 1]]$sqrsum[lb_idxs]
              list_motifs_profile[seq[i], "lb_distances", ] <- lb_profile[lb_idxs]
              list_motifs_profile[seq[i], "index_query", ] <- seq[i]
              list_motifs_profile[seq[i], "indexes_data", ] <- lb_idxs
              list_motifs_profile[seq[i], "dps", ] <- last_product[lb_idxs]
            }
          }

          list_valid_entries <- c(list_valid_entries, indexes_to_update)
          index_valid_entries <- c(index_valid_entries, rep(1, length(indexes_to_update)))
        }

        if (any(!non_valid_idxs)) {
          non_valid_smaller <- min(list_lb_min_non_valid[!non_valid_idxs])
        } else {
          non_valid_smaller <- -1
        }

        if (verbose > 1) {
          pbp$tick(0, tokens = list(what = "Re Check "))
        }

        # Re-check
        m <- (list_motifs_profile[list_valid_entries, "distances", , drop = FALSE][cbind(seq_along(index_valid_entries), 1, index_valid_entries)] < non_valid_smaller)
        motifs_per_size <- sum(m)

        if (motifs_per_size > 0) {
          i <- list_valid_entries[m]
          j <- index_valid_entries[m]
          distances <- list_motifs_profile[i, "distances", , drop = FALSE][cbind(seq_along(j), 1, j)]

          # UPDATE VALMAP_t!
          real_distance <- sqrt(distances)
          normalized_distance <- sqrt(distances) * sqrt(1.0 / window_size)

          n <- (normalized_distance < valmp$matrix_profile[i])

          if (any(n)) {
            entries <- i[n]
            indexes <- list_motifs_profile[entries, "indexes_data", , drop = FALSE][cbind(seq_along(j[n]), 1, j[n])]

            valmp$matrix_profile[entries] <- normalized_distance[n]
            valmp$profile_index[entries] <- indexes
            valmp$length_profile[entries] <- window_size
          }

          r <- (real_distance < valmp$matrix_profile_non_length_normalized[i])

          if (any(r)) {
            entries <- i[r]
            indexes <- list_motifs_profile[entries, "indexes_data", , drop = FALSE][cbind(seq_along(j[r]), 1, j[r])]

            valmp$matrix_profile_non_length_normalized[entries] <- real_distance[r]
            valmp$profile_index_non_length_normalized[entries] <- indexes
            valmp$length_profile_non_length_normalized[entries] <- window_size
          }

          # TODO: check if offset MUST be reduced when re-STOMP

          best_motif <- min(distances)
        }
      }

      if (!is.null(best_motif)) {
        valmp$distances_evolution_motif[offset + 1, ] <- best_motif * sqrt(1.0 / window_size)
      }

      if (max_number_motifs_found < motifs_per_size) {
        max_number_motifs_found <- motifs_per_size
      }
      if ((min_number_motifs_found > motifs_per_size || min_number_motifs_found < 0)) {
        min_number_motifs_found <- motifs_per_size
      }

      if (verbose > 1) {
        pbp$tick(tokens = list(what = "Pruning "))
      }
    }
  }
  tictac <- Sys.time() - tictac

  if (verbose > 1) {
    if (!pbp$finished) {
      pbp$terminate()
    }
  }

  if (verbose > 0) {
    message(sprintf("max_number_motifs_found %s", max_number_motifs_found))
    message(sprintf("min_number_motifs_found %s", min_number_motifs_found))
    message(sprintf("Finished in %.2f %s", tictac, units(tictac)))
  }

  return({
    obj <- list(
      mp = valmp$matrix_profile, pi = valmp$profile_index,
      rmp = NULL, rpi = NULL,
      lmp = NULL, lpi = NULL,
      w = valmp$length_profile,
      ez = ez,
      evolution_motif = valmp$distances_evolution_motif,
      mpnn = valmp$matrix_profile_non_length_normalized,
      pinn = valmp$profile_index_non_length_normalized,
      wnn = valmp$length_profile_non_length_normalized
    )
    class(obj) <- c("Valmod", "MatrixProfile")
    attr(obj, "join") <- join
    obj
  })
}
