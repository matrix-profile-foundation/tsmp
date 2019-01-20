# compute STOMP for lmin and store best 'p' LB DPs /// STOMP_MotifsProfile_VALMOD()
# update VALMP /// matrixProfileVL->distancesEvolutionMotif[offset] = bestMP.distance * ln1SQ;
# for each other window size:
# Compute SubMP
# if best
# update VALMP
# else
# compute STOMP and store best 'p' LB DPs
# update VALMP


# motifIdx is g_orderedidP.bestIdx
# motifDistance is g_orderedidP.bestDist

# if (offset == 0 || motifsPerSize == 0)
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
# realDistance <- sqrt(motifDistance)
# normalizedDistance <- realDistance * sqrt(1.0 / sizeQuery)
# lowerbound <- lbDistance * stdDevQuery^2 / stdDevQueryNext^2
#
# if (normalizedDistance < valmp.matrixProfile || offset == 0)
#   valmp.matrixProfile <- normalizedDistance
#   valmp.indexProfile <- motifIdx # g_orderedidP.bestIdx
#   valmp.lengthProfile <- sizeQuery
#
# if (realDistance < valmp.matrixProfileNonLengthNormalized || offset == 0)
#   valmp.matrixProfileNonLengthNormalized <- realDistance
#   valmp.indexProfileNonLengthNormalized <- motifIdx
#   valmp.lengthProfileNonLengthNormalized <- sizeQuery
#
# stomp_valmod returns bestMP which stores:
# bestMP.distance # motifDistance
# bestMP.index1 # queryIdx
# bestMP.index2 # motifIdx
#
# valmp.distancesEvolutionMotif[offset] <- bestMP.distance * sqrt(1.0 / sizeQuery) # normalizedDistance

valmod <- function(..., window_min, window_max, exclusion_zone = 1 / 2, verbose = 2) {
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
  heap_size <- 150

  if (window_min > query_size / 2) {
    stop("Error: Time series is too short relative to desired window size.", call. = FALSE)
  }
  if (window_min < 4) {
    stop("Error: `window_size` must be at least 4.", call. = FALSE)
  }

  # check skip position
  skip_location <- rep(FALSE, max_profile_size)

  # fixme skip locations
  # for (i in 1:max_profile_size) {
  #   if (any(is.na(data[i:(i + window_min - 1)])) || any(is.infinite(data[i:(i + window_min - 1)]))) {
  #     skip_location[i] <- TRUE
  #   }
  # }

  data[is.na(data)] <- 0
  data[is.infinite(data)] <- 0

  query[is.na(query)] <- 0
  query[is.infinite(query)] <- 0

  valmp <- list(
    matrix_profile = matrix(Inf, max_profile_size, 1),
    profile_index = matrix(-1, max_profile_size, 1),
    length_profile = matrix(-1, max_profile_size, 1),
    distances_evolution_motif = matrix(Inf, range_size, 1),
    matrix_profile_non_length_normalized = matrix(Inf, max_profile_size, 1),
    profile_index_non_length_normalized = matrix(-1, max_profile_size, 1),
    length_profile_non_length_normalized = matrix(-1, max_profile_size, 1)
  )

  list_motifs_profile <- vector(mode = "list", length = max_profile_size)


  for (window_size in seq(window_min, window_max, 1)) {
    # ============== Valmod loop ===============

    message(sprintf("Window size %s", window_size))

    exclusion_zone <- round(window_size * ez + vars()$eps)
    matrix_profile_size <- data_size - window_size + 1
    num_queries <- query_size - window_size + 1
    offset <- window_size - window_min
    motifs_per_size <- 0

    message(sprintf("matrix_profile_size %s", matrix_profile_size))

    if (offset == 0 || motifs_per_size == 0) {
      # ==== STOMP ====

      message(sprintf("Offset %s", offset))

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

      tictac <- Sys.time()

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
      lb_profile <- matrix(0, matrix_profile_size, 1)
      last_product <- matrix(0, matrix_profile_size, 1)
      drop_value <- matrix(0, 1, 1)

      message(sprintf("matrix_profile size %s", length(matrix_profile)))
      message(sprintf("distance_profile size %s", length(distance_profile)))

      for (i in 1:num_queries) {
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

        lb_profile[, 1] <- Re(((last_product / window_size) - (query_mean * data_mean)) / (query_sd * data_sd))
        lb_idx <- (lb_profile > 0)
        lb_profile[lb_idx] <- window_size * (1 - (lb_profile[lb_idx]^2))
        lb_profile[!lb_idx] <- window_size

        drop_value <- query_window[1, 1]

        # apply exclusion zone
        if (exclusion_zone > 0) {
          exc_st <- max(1, i - exclusion_zone)
          exc_ed <- min(matrix_profile_size, i + exclusion_zone)
          distance_profile[exc_st:exc_ed, 1] <- Inf
          lb_profile[exc_st:exc_ed] <- Inf
        }

        distance_profile[data_sd < vars()$eps] <- Inf
        lb_profile[data_sd < vars()$eps] <- Inf
        if (skip_location[i] || any(query_sd[i] < vars()$eps)) {
          distance_profile[] <- Inf
          lb_profile[] <- Inf
        }

        # distance_profile[skip_location] <- Inf
        # lb_profile[skip_location] <- Inf

        lb_idxs <- sort(lb_profile, index.return = TRUE)$ix[1:heap_size]

        list_motifs_profile[[i]] <- list(
          distance = distance_profile[lb_idxs],
          query_sd = query_sd[lb_idxs],
          lb_distance = lb_profile[lb_idxs],
          index1 = i,
          index2 = lb_idxs,
          dp = last_product[lb_idxs]
        )

        distance_profile <- sqrt(distance_profile)

        if (length(args) == 1) {
          # no RMP and LMP for joins
          # left matrix_profile
          ind <- (distance_profile[i:matrix_profile_size] < left_matrix_profile[i:matrix_profile_size])
          ind <- c(rep(FALSE, (i - 1)), ind) # pad left
          left_matrix_profile[ind] <- distance_profile[ind]
          left_profile_index[which(ind)] <- i

          # right matrix_profile
          ind <- (distance_profile[1:i] < right_matrix_profile[1:i])
          ind <- c(ind, rep(FALSE, matrix_profile_size - i)) # pad right
          right_matrix_profile[ind] <- distance_profile[ind]
          right_profile_index[which(ind)] <- i
        }

        ind <- (distance_profile < matrix_profile)
        matrix_profile[ind] <- distance_profile[ind]
        profile_index[which(ind)] <- i

        if (verbose > 1) {
          pb$tick()
        }
      }

      # update valmp
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

      valmp$distances_evolution_motif[offset + 1] <- min(normalized_distance)
    } else {
      #### LB Pruning ####
    }
  }
  tictac <- Sys.time() - tictac

  if (verbose > 0) {
    message(sprintf("Finished in %.2f %s", tictac, units(tictac)))
  }

  return({
    obj <- list(
      #mp = matrix_profile, pi = profile_index,
      mp = valmp$matrix_profile, pi = valmp$profile_index,
      rmp = right_matrix_profile, rpi = right_profile_index,
      lmp = left_matrix_profile, lpi = left_profile_index,
      w = window_size,
      ez = ez,
      valmp = valmp,
      listmotif = list_motifs_profile
    )
    class(obj) <- "MatrixProfile"
    obj
  })
}
