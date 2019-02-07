# compute STOMP for lmin and store best 'p' LB DPs /// STOMP_MotifsProfile_VALMOD()
# update VALMP /// valmp->distancesEvolutionMotif[offset] = bestMP.distance * ln1SQ;
# for each other window size:
# Compute SubMP
# if best
# update VALMP
# else
# compute STOMP and store best 'p' LB DPs
# update VALMP

# motifIdx is g_orderedidP.bestIdx
# motifDistance is g_orderedidP.bestDist

# if (offset == 0 || motifs_per_size == 0)
#    do stomp_valmod
#
# Retrieve the best motif distance and index
# pushes lbDistance, motifDistance, queryIdx, motifIdx, querySize into a list and sort by lbDistance
# lbDistance:
# q <- ((z[i] / querySize) - ((sumQuery / querySize) * (sumData / querySize))) / (stdDevQuery * stdDevData);
# if (q > 0)
#   lbDistance <- querySize * (1 - (q * q))
# else
#   lbDistance <- querySize
#
# real_distance <- sqrt(motifDistance)
# normalized_distance <- real_distance * sqrt(1.0 / sizeQuery)
# lowerbound <- lbDistance * stdDevQuery^2 / stdDevQueryNext^2
#
# if (normalized_distance < valmp.matrix_profile || offset == 0)
#   valmp.matrix_profile <- normalized_distance
#   valmp.profile_index <- motifIdx # g_orderedidP.bestIdx
#   valmp.length_profile <- sizeQuery
#
# if (real_distance < valmp.matrix_profile_non_length_normalized || offset == 0)
#   valmp.matrix_profile_non_length_normalized <- real_distance
#   valmp.profile_index_non_length_normalized <- motifIdx
#   valmp.length_profile_non_length_normalized <- sizeQuery
#
# stomp_valmod returns bestMP which stores:
# bestMP.distance # motifDistance
# bestMP.index_query # queryIdx
# bestMP.indexes_data # motifIdx
#
# valmp.distancesEvolutionMotif[offset] <- bestMP.distance * sqrt(1.0 / sizeQuery) # normalized_distance
#' @export

valmod <- function(..., window_min, window_max, heap_size = 50, exclusion_zone = 1 / 2, verbose = 2, lb = FALSE) {
  args <- list(...)
  data <- args[[1]]
  if (length(args) > 1) {
    query <- args[[2]]
    exclusion_zone <- 0 # don't use exclusion zone for joins
  } else {
    query <- data
  }

  # transform data into matrix
  if (is.vector(data)) {
    data <- as.matrix(data)
  }
  else if (is.matrix(data)) {
    if (ncol(data) > nrow(data)) {
      data <- t(data)
    }
  } else {
    stop("Error: Unknown type of data. Must be: a column matrix or a vector.", call. = FALSE)
  }

  if (is.vector(query)) {
    query <- as.matrix(query)
  } else if (is.matrix(query)) {
    if (ncol(query) > nrow(query)) {
      query <- t(query)
    }
  } else {
    stop("Error: Unknown type of query. Must be: a column matrix or a vector.", call. = FALSE)
  }

  ez <- exclusion_zone # store original
  data_size <- nrow(data)
  query_size <- nrow(query)
  range_size <- window_max - window_min + 1
  max_profile_size <- data_size - window_min + 1

  if (window_min > query_size / 2) {
    stop("Error: Time series is too short relative to desired window size.", call. = FALSE)
  }
  if (window_min < 4) {
    stop("Error: `window_size` must be at least 4.", call. = FALSE)
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


  matrix_profiles_elements_per_size <- 1
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

  for (window_size in seq(window_min, window_max, 1)) {
    # ============== Valmod loop ===============
    message(sprintf("Window size %s", window_size))

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
      message("STOMP")
      tictac_stomp <- Sys.time()

      if (verbose > 1) {
        pb <- progress::progress_bar$new(
          format = "STOMP [:bar] :percent at :tick_rate it/s, elapsed: :elapsed, eta: :eta",
          clear = FALSE, total = num_queries, width = 80
        )
      }

      if (verbose > 2) {
        on.exit(beep(sounds[[1]]), TRUE)
      }

      data_fft <- matrix(0, (window_size + data_size), 1)
      data_mean <- matrix(0, matrix_profile_size, 1)
      data_sd <- matrix(0, matrix_profile_size, 1)
      first_product <- matrix(0, num_queries, 1)

      # forward
      nnpre <- mass_pre(data, data_size, query, query_size, window_size = window_size)
      data_fft <- nnpre$data_fft
      data_mean <- nnpre$data_mean
      data_sd <- nnpre$data_sd
      query_mean <- nnpre$query_mean
      query_sd <- nnpre$query_sd

      # reverse
      # This is needed to handle with the join similarity.
      rnnpre <- mass_pre(query, query_size, data, data_size, window_size = window_size)
      rdata_fft <- rnnpre$data_fft
      rdata_mean <- rnnpre$data_mean
      rdata_sd <- rnnpre$data_sd
      rquery_mean <- rnnpre$query_mean
      rquery_sd <- rnnpre$query_sd

      rnn <- mass(
        rdata_fft, data[1:window_size], query_size, window_size, rdata_mean,
        rdata_sd, rquery_mean[1], rquery_sd[1]
      )

      first_product[, 1] <- rnn$last_product

      matrix_profile <- matrix(Inf, matrix_profile_size, 1)
      profile_index <- matrix(-1, matrix_profile_size, 1)
      if (length(args) > 1) {
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
          nn <- mass(
            data_fft, query_window, data_size, window_size, data_mean, data_sd,
            query_mean[i], query_sd[i]
          )
          distance_profile[, 1] <- nn$distance_profile
          last_product[, 1] <- nn$last_product
        } else {
          last_product[2:(data_size - window_size + 1), 1] <- last_product[1:(data_size - window_size), 1] -
            data[1:(data_size - window_size), 1] * drop_value +
            data[(window_size + 1):data_size, 1] * query_window[window_size, 1]

          last_product[1, 1] <- first_product[i, 1]
          distance_profile <- 2 * (window_size - (last_product - window_size * data_mean * query_mean[i]) /
            (data_sd * query_sd[i]))
        }

        # distance_profile <- Re(sqrt(distance_profile))
        distance_profile <- Re(distance_profile)
        distance_profile[distance_profile < 0] <- 0

        new_lb_profile <- Re(((last_product / window_size) - (query_mean[i] * data_mean)) / (query_sd[i] * data_sd))
        new_lb_profile[new_lb_profile < 0] <- 0
        new_lb_profile_len <- length(new_lb_profile)

        if (new_lb_profile_len != matrix_profile_size) {
          stop("new_lb_profile_len != matrix_profile_size")
        }

        lb_profile <- new_lb_profile
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

        distance_profile[data_sd < vars()$eps] <- Inf
        lb_profile[data_sd < vars()$eps] <- Inf
        if (skip_location[i] || any(query_sd[i] < vars()$eps)) {
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
        list_motifs_profile[i, "query_sd", ] <- query_sd[i]
        list_motifs_profile[i, "sum_query", ] <- query_stats[[offset + 1]]$sum[i]
        list_motifs_profile[i, "sum_data", ] <- data_stats[[offset + 1]]$sum[lb_idxs]
        list_motifs_profile[i, "sqrsum_query", ] <- query_stats[[offset + 1]]$sqrsum[i]
        list_motifs_profile[i, "sqrsum_data", ] <- data_stats[[offset + 1]]$sqrsum[lb_idxs]
        list_motifs_profile[i, "lb_distances", ] <- lb_profile[lb_idxs]
        list_motifs_profile[i, "index_query", ] <- i
        list_motifs_profile[i, "indexes_data", ] <- lb_idxs
        list_motifs_profile[i, "dps", ] <- Re(last_product[lb_idxs])

        distance_profile <- sqrt(distance_profile)

        # if (length(args) == 1) {
        #   # no RMP and LMP for joins
        #   # left matrix_profile
        #   ind <- (distance_profile[i:matrix_profile_size] < left_matrix_profile[i:matrix_profile_size])
        #   ind <- c(rep(FALSE, (i - 1)), ind) # pad left
        #   left_matrix_profile[ind] <- distance_profile[ind]
        #   left_profile_index[which(ind)] <- i
        #
        #   # right matrix_profile
        #   ind <- (distance_profile[1:i] < right_matrix_profile[1:i])
        #   ind <- c(ind, rep(FALSE, matrix_profile_size - i)) # pad right
        #   right_matrix_profile[ind] <- distance_profile[ind]
        #   right_profile_index[which(ind)] <- i
        # }

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
      message(sprintf("distances_evolution_motif %s", valmp$distances_evolution_motif[offset + 1, ]))

      motifs_per_size <- -1

      tictac_stomp <- Sys.time() - tictac_stomp

      if (verbose > 0) {
        message(sprintf("Finished STOMP in %.2f %s", tictac_stomp, units(tictac_stomp)))
      }
    } else {
      #### LB Pruning ####
      message("LB Pruning")
      tictac_lb <- Sys.time()

      if (offset == 0) {
        message("OFFSET ZERO?")
      }
      # TODO: offset never gets zero in this block, why check offset == 0 ?

      matrix_profiles_elements_per_size <- 1
      motifs_per_size <- 0

      list_valid_entries <- vector(mode = "numeric", length = matrix_profile_size)
      index_valid_entries <- vector(mode = "numeric", length = matrix_profile_size)
      non_valid_smaller <- -1

      max_entry_valid <- NULL
      min_abs_true_dist <- -1

      non_valid_lb <- 1

      for (i in seq_len(matrix_profile_size)) {
        # --- First Outer Loop ----

        max_lb <- list_motifs_profile[i, "lb_distances", heap_size]
        max_query_sd <- list_motifs_profile[i, "query_sd", heap_size]
        curr_query_sd <- query_stats[[offset + 1]]$sd[i] # std(query[i:(i + window_size - 1)])
        lower_bound <- max_lb * max_query_sd^2 / curr_query_sd^2

        min_entry_idx <- NULL
        min_entry_lb <- -1

        # message(sprintf("heap_size %s", heap_size))
        new_index_query <- list_motifs_profile[i, "index_query", 1] + window_size - 1

        if (new_index_query > query_size) {
          next # not valid offset
        }

        list_motifs_profile[i, "sum_query", ] <- list_motifs_profile[i, "sum_query", ] + query[new_index_query]
        list_motifs_profile[i, "sqrsum_query", ] <- list_motifs_profile[i, "sqrsum_query", ] + (query[new_index_query] * query[new_index_query])
        query_mean_v <- list_motifs_profile[i, "sum_query", ] / window_size
        query_sd_v <- (list_motifs_profile[i, "sqrsum_query", ] / window_size) - (query_mean_v * query_mean_v)

        query_sd_v[query_sd_v < 0] <- 0
        query_sd_v <- sqrt(query_sd_v)

        # --- First Inner Loop ----
        j_v <- seq(heap_size, 2, by = -1)

        # apply exclusion zone
        ez_v <- (list_motifs_profile[i, "indexes_data", j_v] <= max(1, list_motifs_profile[i, "indexes_data", j_v] - exclusion_zone) |
          list_motifs_profile[i, "indexes_data", j_v] >= min(matrix_profile_size, list_motifs_profile[i, "indexes_data", j_v] + exclusion_zone)) |
          list_motifs_profile[i, "indexes_data", j_v] + window_size - 1 <= data_size
        j_v <- j_v[ez_v]

        new_index_data_v <- list_motifs_profile[i, "indexes_data", j_v] + window_size - 1

        # --- Compute true distance entry heapentry ----
        list_motifs_profile[i, "dps", j_v] <- list_motifs_profile[i, "dps", j_v] + (query[new_index_query] * data[new_index_data_v])
        list_motifs_profile[i, "sum_data", j_v] <- list_motifs_profile[i, "sum_data", j_v] + data[new_index_data_v]
        list_motifs_profile[i, "sqrsum_data", j_v] <- list_motifs_profile[i, "sqrsum_data", j_v] + (data[new_index_data_v] * data[new_index_data_v])

        data_mean_v <- list_motifs_profile[i, "sum_data", j_v] / window_size
        data_sd_v <- (list_motifs_profile[i, "sqrsum_data", j_v] / window_size) - (data_mean_v * data_mean_v)

        data_sd_v[data_sd_v < 0] <- 0
        data_sd_v <- sqrt(data_sd_v)

        dist_v <- ((2 * window_size) * (1 - ((list_motifs_profile[i, "dps", j_v] - (window_size * query_mean_v[j_v] * data_mean_v)) /
          (window_size * query_sd_v[j_v] * data_sd_v))))

        dist_v <- Re(dist_v)
        dist_v[dist_v < 0 ] <- 0

        list_motifs_profile[i, "distances", j_v] <- dist_v

        # ---- Define min_entry ----
        min_entry_idx <- j_v[which.min(list_motifs_profile[i, "distances", j_v])]
        min_entry_lb <- max(list_motifs_profile[i, "lb_distances", j_v] * list_motifs_profile[i, "query_sd", j_v]^2 / curr_query_sd^2)

        # if (is.null(min_entry_idx)) {
        #   min_entry_idx <- j_v
        #   min_entry_lb <- list_motifs_profile[i, "lb_distances", j_v] * list_motifs_profile[i, "query_sd", j_v]^2 / curr_query_sd^2
        # }
        # else {
        #   if (list_motifs_profile[i, "distances", min_entry_idx] > dist_v) {
        #     min_entry_idx <- j_v
        #   }

        #   lb <- list_motifs_profile[i, "lb_distances", j_v] * list_motifs_profile[i, "query_sd", j_v]^2 / curr_query_sd^2
        #   if (min_entry_lb < lb) {
        #     min_entry_lb <- lb
        #   }
        # }

        if (!is.null(min_entry_idx) && list_motifs_profile[i, "distances", min_entry_idx] < lower_bound) {
          #### Pruning is effective global min found ####
          # update the min OOOOONLYYYY IF it is a correct min of the distance profile

          if (min_abs_true_dist == -1) {
            min_abs_true_dist <- list_motifs_profile[i, "distances", min_entry_idx]
          } else if (list_motifs_profile[i, "distances", min_entry_idx] < min_abs_true_dist) {
            min_abs_true_dist <- list_motifs_profile[i, "distances", min_entry_idx]
          }

          list_valid_entries[matrix_profiles_elements_per_size] <- i
          index_valid_entries[matrix_profiles_elements_per_size] <- min_entry_idx

          # new element of the matrix profile
          matrix_profiles_elements_per_size <- matrix_profiles_elements_per_size + 1
          if (is.null(max_entry_valid)) {
            max_entry_valid <- list_motifs_profile[i, "distances", min_entry_idx]
          } else if (max_entry_valid < list_motifs_profile[i, "distances", min_entry_idx]) {
            max_entry_valid <- list_motifs_profile[i, "distances", min_entry_idx]
          }
        } else if (min_entry_lb >= 0) {
          #### The distance profile has elements but nothing to say, store the min max LB ####
          list_lb_min_non_valid[non_valid_lb] <- min_entry_lb
          index_lb_min_non_valid[non_valid_lb] <- i

          # store the non valid entry with the smallest LB.
          # (I am sure that this guy correctly lowerbounds the MP  )

          if (non_valid_smaller == -1) {
            non_valid_smaller <- min_entry_lb
          } else if (non_valid_smaller > min_entry_lb) {
            non_valid_smaller <- min_entry_lb
          }

          non_valid_lb <- non_valid_lb + 1
        } else {
          #### All trival matches recompute the DP ####
          list_lb_min_non_valid[non_valid_lb] <- -1
          index_lb_min_non_valid[non_valid_lb] <- i
        }
      }

      # check to see how many valid motifs (smallest value of the Matrix Profile) I have among
      # the valid entries.

      best_motif <- NULL

      message(sprintf("matrix_profiles_elements_per_size %s", matrix_profiles_elements_per_size))

      for (m in seq_len(matrix_profiles_elements_per_size - 1)) {
        # --- Second Outer Loop ----

        i <- list_valid_entries[m]
        j <- index_valid_entries[m]

        if (list_motifs_profile[i, "distances", j] < non_valid_smaller) {
          # UPDATE VALMAP_t!
          real_distance <- sqrt(list_motifs_profile[i, "distances", j])
          normalized_distance <- sqrt(list_motifs_profile[i, "distances", j]) * sqrt(1.0 / window_size)

          index_update_vm <- list_motifs_profile[i, "index_query", j]

          if (normalized_distance < valmp$matrix_profile[index_update_vm] || offset == 0) {
            valmp$matrix_profile[index_update_vm] <- normalized_distance
            valmp$profile_index[index_update_vm] <- list_motifs_profile[i, "indexes_data", j]
            valmp$length_profile[index_update_vm] <- window_size
          }

          if (real_distance < valmp$matrix_profile_non_length_normalized[index_update_vm] ||
            offset == 0) {
            valmp$matrix_profile_non_length_normalized[index_update_vm] <- real_distance
            valmp$profile_index_non_length_normalized[index_update_vm] <- list_motifs_profile[i, "indexes_data", j]
            valmp$length_profile_non_length_normalized[index_update_vm] <- window_size
          }

          # TODO: check why distancesEvolutionMotif is not updated
          # TODO: check if offset MUST be reduced when re-STOMP

          if (is.null(best_motif)) {
            best_motif <- list_motifs_profile[i, "distances", j]
          } else if (best_motif > list_motifs_profile[i, "distances", j]) {
            best_motif <- list_motifs_profile[i, "distances", j]
          }

          motifs_per_size <- motifs_per_size + 1
        }
      }

      if (is.null(best_motif)) {
        message("best_motif NULL")

        # matrix profile may not be containing the minimum
        # Here I do not want to call STOMP just update the distance profiles which have the max
        # LB smaller than the minimum motif pair

        min_lb <- -1

        message(sprintf("non_valid_lb %s", non_valid_lb))

        for (z in seq_len(non_valid_lb - 1)) {
          if (list_lb_min_non_valid[z] <= min_abs_true_dist) {
            index_to_update <- index_lb_min_non_valid[z]

            if (index_to_update > matrix_profile_size) {
              stop("index_to_update > matrix_profile_size")
            }

            pre <- mass_pre(data, data_size, query, query_size, window_size = window_size)

            nn <- mass(
              pre$data_fft, query[index_to_update:(index_to_update + window_size - 1)], data_size,
              window_size, pre$data_mean, pre$data_sd, pre$query_mean[index_to_update],
              pre$query_sd[index_to_update]
            )

            distance_profile <- Re(nn$distance_profile)
            distance_profile[distance_profile < 0] <- 0

            new_lb_profile <- Re(((nn$last_product / window_size) - (pre$query_mean[index_to_update] *
              pre$data_mean)) / (pre$query_sd[index_to_update] * pre$data_sd))
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
              exc_st <- max(1, index_to_update - exclusion_zone)
              exc_ed <- min(matrix_profile_size, index_to_update + exclusion_zone)

              distance_profile[exc_st:exc_ed] <- Inf
              lb_profile[exc_st:exc_ed] <- Inf
            }

            distance_profile[pre$data_sd < vars()$eps] <- Inf
            lb_profile[pre$data_sd < vars()$eps] <- Inf
            if (skip_location[index_to_update] || any(pre$query_sd[index_to_update] < vars()$eps)) {
              distance_profile[] <- Inf
              lb_profile[] <- Inf
            }

            distance_profile[skip_location] <- Inf
            lb_profile[skip_location] <- Inf

            lb_idxs <- sort(lb_profile, index.return = TRUE)$ix[1:heap_size]

            if (any(lb_idxs > matrix_profile_size)) {
              stop("lb_idxs > matrix_profile_size")
            }

            list_motifs_profile[index_to_update, "distances", ] <- distance_profile[lb_idxs]
            list_motifs_profile[index_to_update, "query_sd", ] <- query_sd_v
            list_motifs_profile[index_to_update, "sum_query", ] <- query_stats[[offset + 1]]$sum[index_to_update]
            list_motifs_profile[index_to_update, "sum_data", ] <- data_stats[[offset + 1]]$sum[lb_idxs]
            list_motifs_profile[index_to_update, "sqrsum_query", ] <- query_stats[[offset + 1]]$sqrsum[index_to_update]
            list_motifs_profile[index_to_update, "sqrsum_data", ] <- data_stats[[offset + 1]]$sqrsum[lb_idxs]
            list_motifs_profile[index_to_update, "lb_distances", ] <- lb_profile[lb_idxs]
            list_motifs_profile[index_to_update, "index_query", ] <- index_to_update
            list_motifs_profile[index_to_update, "indexes_data", ] <- lb_idxs
            list_motifs_profile[index_to_update, "dps", ] <- Re(last_product[lb_idxs])

            list_valid_entries[matrix_profiles_elements_per_size] <- index_to_update
            index_valid_entries[matrix_profiles_elements_per_size] <- 1

            matrix_profiles_elements_per_size <- matrix_profiles_elements_per_size + 1
          } else {
            if (min_lb < 0) {
              min_lb <- list_lb_min_non_valid[z]
            } else if (min_lb > list_lb_min_non_valid[z]) {
              min_lb <- list_lb_min_non_valid[z]
            }
          }
        }

        non_valid_smaller <- min_lb

        message(sprintf("matrix_profiles_elements_per_size %s", matrix_profiles_elements_per_size))

        # Re-check
        for (m in seq_len(matrix_profiles_elements_per_size - 1)) {
          i <- list_valid_entries[m]
          j <- index_valid_entries[m]

          if (list_motifs_profile[i, "distances", j] < non_valid_smaller) {
            # UPDATE VALMAP_t!
            real_distance <- sqrt(list_motifs_profile[i, "distances", j])
            normalized_distance <- sqrt(list_motifs_profile[i, "distances", j]) * sqrt(1.0 / window_size)

            index_update_vm <- list_motifs_profile[i, "index_query", j]

            if (index_update_vm > matrix_profile_size) {
              stop("index_update_vm > matrix_profile_size")
            }

            if (normalized_distance < valmp$matrix_profile[index_update_vm] || offset == 0) {
              valmp$matrix_profile[index_update_vm] <- normalized_distance
              valmp$profile_index[index_update_vm] <- list_motifs_profile[i, "indexes_data", j]
              valmp$length_profile[index_update_vm] <- window_size
            }
            if (real_distance < valmp$matrix_profile[index_update_vm] || offset == 0) {
              valmp$matrix_profile_non_length_normalized[index_update_vm] <- real_distance
              valmp$profile_index_non_length_normalized[index_update_vm] <- list_motifs_profile[i, "indexes_data", j]
              valmp$length_profile_non_length_normalized[index_update_vm] <- window_size
            }

            if (is.null(best_motif)) {
              best_motif <- list_motifs_profile[i, "distances", j]
            } else {
              if (best_motif > list_motifs_profile[i, "distances", j]) {
                best_motif <- list_motifs_profile[i, "distances", j]
              }
            }
            motifs_per_size <- motifs_per_size + 1
          }
        }
      }

      if (!is.null(best_motif)) {
        valmp$distances_evolution_motif[offset + 1, ] <- best_motif * sqrt(1.0 / window_size)
        message(sprintf("distances_evolution_motif %s", valmp$distances_evolution_motif[offset + 1, ]))
      }

      if (max_number_motifs_found < motifs_per_size) {
        max_number_motifs_found <- motifs_per_size
      }
      if ((min_number_motifs_found > motifs_per_size || min_number_motifs_found < 0)) {
        min_number_motifs_found <- motifs_per_size
      }

      tictac_lb <- Sys.time() - tictac_lb

      if (verbose > 0) {
        message(sprintf("Finished LB in %.2f %s", tictac_lb, units(tictac_lb)))
      }
    }
  }
  tictac <- Sys.time() - tictac

  if (verbose > 0) {
    message(sprintf("Finished in %.2f %s", tictac, units(tictac)))
  }

  return({
    obj <- list(
      mp = valmp$matrix_profile, pi = valmp$profile_index,
      rmp = right_matrix_profile, rpi = right_profile_index,
      lmp = left_matrix_profile, lpi = left_profile_index,
      w = valmp$length_profile,
      valmp = valmp,
      ez = ez
    )
    class(obj) <- "matrix_profile"
    obj
  })
}
