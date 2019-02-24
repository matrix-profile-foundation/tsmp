#' Search for Discord
#'
#' @param .mp a TSMP object of class `MatrixProfile`
#' @param ... further arguments to be passed to class specific function.
#' @name find_discord
#' @export

find_discord <- function(.mp, ...) {
  UseMethod("find_discord", .mp)
}

#' @param data the data used to build the Matrix Profile, if not embedded.
#' @param n_discords an `int`. Number of discords to find. (Default is `1`).
#' @param n_neighbors an `int`. Number of neighbors to find. (Default is `3`).
#' @param radius an `int`. Set a threshold to exclude matching neighbors with distance > current
#' discord distance * `radius`. (Default is `3`).
#' @param exclusion_zone if a `number` will be used instead of embedded value. (Default is `NULL`).
#' @name find_discord
#' @export
#' @return For class `MatrixProfile`, returns the input `.mp` object with a new name `discord`. It contains: `discord_idx`, a `vector`
#' of discords founded.
#' @examples
#' # Single dimension data
#' w <- 50
#' data <- mp_gait_data
#' mp <- tsmp(data, window_size = w, exclusion_zone = 1 / 4, verbose = 0)
#' mp <- find_discord(mp)
find_discord.MatrixProfile <- function(.mp, data, n_discords = 1, n_neighbors = 3, radius = 3, exclusion_zone = NULL, ...) {
  if (!("MatrixProfile" %in% class(.mp))) {
    stop("Error: First argument must be an object of class `MatrixProfile`.")
  }

  if ("Valmod" %in% class(.mp)) {
    stop("Error: Function not implemented for objects of class `Valmod`.")
  }

  if (missing(data) && !is.null(.mp$data)) {
    data <- .mp$data[[1]]
  }

  # transform data list into matrix
  if (is.matrix(data) || is.data.frame(data)) {
    if (is.data.frame(data)) {
      data <- as.matrix(data)
    } # just to be uniform
    if (ncol(data) > nrow(data)) {
      data <- t(data)
    }
    data_len <- nrow(data)
    data_dim <- ncol(data)
  } else if (is.list(data)) {
    data_len <- length(data[[1]])
    data_dim <- length(data)

    for (i in 1:data_dim) {
      len <- length(data[[i]])
      # Fix TS size with NaN
      if (len < data_len) {
        data[[i]] <- c(data[[i]], rep(NA, data_len - len))
      }
    }
    # transform data into matrix (each column is a TS)
    data <- sapply(data, cbind)
  } else if (is.vector(data)) {
    data_len <- length(data)
    data_dim <- 1
    # transform data into 1-col matrix
    data <- as.matrix(data) # just to be uniform
  } else {
    stop("Error: `data` must be `matrix`, `data.frame`, `vector` or `list`.")
  }

  if (is.null(exclusion_zone)) {
    exclusion_zone <- .mp$ez
  }

  matrix_profile <- .mp$mp # keep mp intact
  exclusion_zone <- round(.mp$w * exclusion_zone + vars()$eps)
  data_size <- nrow(data)
  matrix_profile_size <- length(matrix_profile)
  discord_idxs <- list(discords = vector(mode = "numeric"), neighbors = list(NULL))

  nn_pre <- mass_pre(data, data_size, window_size = .mp$w)

  for (i in seq_len(n_discords)) {
    discord_idx <- which.max(matrix_profile)
    discord_distance <- matrix_profile[discord_idx]
    discord_idxs[[1]][i] <- discord_idx

    # query using the discord to find its neighbors
    query <- data[discord_idx:(discord_idx + .mp$w - 1)]

    distance_profile <- mass(
      nn_pre$data_fft, query, data_size, .mp$w, nn_pre$data_mean, nn_pre$data_sd,
      nn_pre$data_mean[discord_idx], nn_pre$data_sd[discord_idx]
    )

    distance_profile <- Re(distance_profile$distance_profile)
    distance_profile[distance_profile > (discord_distance * radius)^2] <- Inf
    discord_zone_start <- max(1, discord_idx - exclusion_zone)
    discord_zone_end <- min(matrix_profile_size, discord_idx + exclusion_zone)
    distance_profile[discord_zone_start:discord_zone_end] <- Inf
    st <- sort(distance_profile, index.return = TRUE)
    distance_order <- st$x
    distance_idx_order <- st$ix

    discord_neighbor <- vector(mode = "numeric")

    for (j in seq_len(n_neighbors)) {
      if (is.infinite(distance_order[1]) || length(distance_order) < j) {
        break
      }
      discord_neighbor[j] <- distance_idx_order[1]
      distance_order <- distance_order[2:length(distance_order)]
      distance_idx_order <- distance_idx_order[2:length(distance_idx_order)]
      distance_order <- distance_order[!(abs(distance_idx_order - discord_neighbor[j]) < exclusion_zone)]
      distance_idx_order <- distance_idx_order[!(abs(distance_idx_order - discord_neighbor[j]) < exclusion_zone)]
    }

    discord_neighbor <- discord_neighbor[discord_neighbor != 0]
    discord_idxs[[2]][[i]] <- discord_neighbor

    remove_idx <- c(discord_idxs[[1]][i], discord_idxs[[2]][[i]])

    for (j in seq_len(length(remove_idx))) {
      remove_zone_start <- max(1, remove_idx[j] - exclusion_zone)
      remove_zone_end <- min(matrix_profile_size, remove_idx[j] + exclusion_zone)
      matrix_profile[remove_zone_start:remove_zone_end] <- -Inf
    }
  }

  .mp$discord <- list(discord_idx = discord_idxs[[1]], discord_neighbor = discord_idxs[[2]])
  class(.mp) <- update_class(class(.mp), "Discord")
  return(.mp)
}
