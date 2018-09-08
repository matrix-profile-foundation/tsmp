
# Math functions ------------------------------------------------------------------------------

#' Fast implementation of moving standard deviation using filter
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#' @param window_size moving sd window size
#'
#' @return Returns a `vector` with the moving standard deviation
#' @export
#'
#' @examples
#' data_sd <- fast_movsd(mp_toy_data$data[,1], mp_toy_data$sub_len)

fast_movsd <- function(data, window_size) {

  # length of the time series
  data_size <- length(data)

  if (window_size < 2) {
    stop("Error: 'window_size' must be at least 2.")
  }

  if (data_size < window_size) {
    stop("Error: 'window_size' is too large for this series.")
  }

  # Improve the numerical analysis by subtracting off the series mean
  # this has no effect on the standard deviation.
  data <- data - mean(data)

  # scale the data to have unit variance too. will put that
  # scale factor back into the result at the end
  data_sd <- std(data)
  data <- data / data_sd

  # we will need the squared elements
  data_sqr <- data^2

  b <- matrix(1, 1, window_size)
  s <- sqrt((stats::filter(data_sqr, b, sides = 1) - (stats::filter(data, b, sides = 1)^2) * (1 / window_size)) / (window_size - 1))

  # restore the scale factor that was used before to normalize the data
  s <- s * data_sd
  s <- Re(s)
  s <- s * sqrt((window_size - 1) / window_size)

  return(s[!is.na(s)])
}

#' Fast implementation of moving average using filter
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#' @param window_size moving average window size
#'
#' @return Returns a `vector` with the moving average
#' @export
#' @examples
#' data_avg <- fast_movavg(mp_toy_data$data[,1], mp_toy_data$sub_len)

fast_movavg <- function(data, window_size) {
  data_mean <- stats::filter(data, rep(1 / window_size, window_size), sides = 2)
  return(data_mean[!is.na(data_mean)])
}

#' Population SD, as R always calculate with n-1 (sample), here we fix it
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#'
#' @return Returns the corrected standard deviation from sample to population
#' @keywords internal
#' @noRd
#'
std <- function(data) {
  sdx <- stats::sd(data)

  if (sdx == 0) {
    return(sdx)
  }

  return(sqrt((length(data) - 1) / length(data)) * sdx)
}

#' Normalizes data for mean Zero and Standard Deviation One
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#'
#' @return Returns the normalized data
#' @keywords internal
#' @noRd
#'
znorm <- function(data) {
  data_mean <- mean(data)
  data_dev <- std(data)

  if (is.nan(data_dev) || data_dev <= 0.01) {
    return(data - data_mean)
  }
  else {
    (data - data_mean) / (data_dev)
  }
}

#' Distance between two matrices
#'
#' Computes the Euclidean distance between rows of two matrices.
#'
#' @param x a `matrix`.
#' @param y a `matrix`.
#'
#' @return Returns a `matrix` of size m x n if x is of size m x k and y is of size n x k.
#' @keywords internal
#' @noRd

diff2 <- function(x, y) {
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("Error: `x` and `y` must be numeric vectors or matrices.")
  }
  if (is.vector(x)) {
    dim(x) <- c(1, length(x))
  }
  if (is.vector(y)) {
    dim(y) <- c(1, length(y))
  }
  if (ncol(x) != ncol(y)) {
    stop("Error: `x` and `y` must have the same number of columns.")
  }
  m <- nrow(x)
  n <- nrow(y)
  xy <- x %*% t(y)
  xx <- matrix(rep(apply(x * x, 1, sum), n), m, n, byrow = FALSE)
  yy <- matrix(rep(apply(y * y, 1, sum), m), m, n, byrow = TRUE)
  sqrt(pmax(xx + yy - 2 * xy, 0))
}

# Salient Aux functions --------------------------------------------------------------------------

#' Reduced description length
#'
#' @param x the difference between two time series (reference and candidate for compression)
#' @param mismatch_bit sum of n_bits and log2(window_size)
#'
#' @return Returns the bit_size cost of compressing the time series
#' @keywords internal
#' @noRd

get_bitsize <- function(x, mismatch_bit) {
  bit_size <- sum(x != 0) * mismatch_bit

  return(bit_size)
}

#' Precompute the max and min value for the discrete normalization
#'
#' @param data input time series
#' @param window_size sliding window size
#'
#' @return Returns a list with the max and min value
#' @keywords internal
#' @noRd
#'
discrete_norm_pre <- function(data, window_size = 1) {
  if (is.vector(data)) {
    data <- as.matrix(data)
  }

  if (ncol(data) > 1) {
    len <- ncol(data)
  } else {
    len <- nrow(data) - window_size + 1
  }

  max <- -Inf
  min <- Inf
  for (i in 1:len) {
    if (ncol(data) > 1) {
      window <- data[, i]
    } else {
      window <- data[i:(i + window_size - 1), ]
    }
    window_mean <- mean(window)
    window_sd <- std(window)
    if (window_sd == 0) {
      window <- (window - window_mean)
    } else {
      window <- (window - window_mean) / window_sd
    }

    if (max(window) > max) {
      max <- max(window)
    }
    if (min(window) < min) {
      min <- min(window)
    }
  }
  return(list(max = max, min = min))
}


#' Discrete normalization
#'
#' @param data Input time series.
#' @param n_bits Number of bits for MDL discretization.
#' @param max Precomputed max from `discrete_norm_pre`.
#' @param min Precomputed min from `discrete_norm_pre`.
#'
#' @return Returns the data after discrete normalization.
#' @keywords internal
#' @noRd

discrete_norm <- function(data, n_bits, max, min) {
  # normalize magnitude
  data_mean <- mean(data)
  data_sd <- std(data)

  if (data_sd == 0) {
    data <- (data - data_mean)
  } else {
    data <- (data - data_mean) / data_sd
  }

  data <- (data - min) / (max - min)

  # discretization
  data <- round(data * (2^n_bits - 1) + vars()$eps) + 1

  return(data)
}

# AV Aux functions --------------------------------------------------------------------------------

#' Counts number of zero-crossings
#'
#' Count the number of zero-crossings from the supplied time-domain input vector. A simple method is
#' applied here that can be easily ported to a real-time system that would minimize the number of
#' if-else conditionals.
#'
#' @param data is the input time-domain signal (one dimensional).
#'
#' @return Returns the amount of zero-crossings in the input signal.
#' @author sparafucile17 06/27/04
#' @references <https://www.dsprelated.com/showcode/179.php>
#' @keywords internal
#' @noRd
#'
zero_crossings <- function(data) {
  # initial value
  count <- 0

  data <- as.matrix(data)

  # error checks
  if (length(data) == 1) {
    stop("Error: Input signal must have more than one element.")
  }

  if ((ncol(data) != 1) && (nrow(data) != 1)) {
    stop("Error: Input must be one-dimensional.")
  }

  # force signal to be a vector oriented in the same direction
  data <- as.vector(data)

  num_samples <- length(data)

  for (i in 2:num_samples) {
    # Any time you multiply to adjacent values that have a sign difference
    # the result will always be negative.  When the signs are identical,
    # the product will always be positive.
    if ((data[i] * data[i - 1]) < 0) {
      count <- count + 1
    }
  }

  return(count)
}

#' Normalizes data between Zero and One
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#'
#' @return Returns the normalized data.
#' @keywords internal
#' @noRd
#'
zero_one_norm <- function(data) {
  data <- round(data, 10)

  data <- data - min(data[!is.infinite(data) & !is.na(data)])
  data <- data / max(data[!is.infinite(data) & !is.na(data)])

  return(data)
}

#' Computes the complexity index of the data
#'
#' @param data a `vector` or a column `matrix` of `numeric`.
#'
#' @return Returns the complexity index of the data provided (normally a subset)
#' @keywords internal
#' @noRd
#'
complexity <- function(data) {
  return(sqrt(sum(diff(data)^2)))
}

# Find Motif Ung Aux Functions -------------------------------------------------------------

#' Computes the bits saved using basic compactation algorithm
#'
#' @param motif_1 reference motif
#' @param motif_2 difference motif
#' @param n_dim data dimensions
#' @param n_bit bits for discretization
#'
#' @return Returns the amount of bits saved
#' @keywords internal
#' @noRd
#'
get_bit_save <- function(motif_1, motif_2, n_dim, n_bit) {
  if (is.vector(motif_1)) {
    motif_1 <- as.matrix(motif_1)
  }

  if (is.vector(motif_2)) {
    motif_2 <- as.matrix(motif_2)
  }

  tot_dim <- dim(motif_1)[2]
  window_size <- dim(motif_1)[1]
  split_pt <- get_desc_split_pt(n_bit)
  disc_1 <- discretization(motif_1, split_pt)
  disc_2 <- discretization(motif_2, split_pt)

  dim_id <- sort(apply(abs(disc_1 - disc_2), 2, sum), index.return = TRUE)$ix
  dim_id <- dim_id[1:n_dim]
  motif_diff <- disc_1[, dim_id] - disc_2[, dim_id]
  n_val <- length(unique(as.vector(motif_diff)))

  bit_sz <- n_bit * (tot_dim * window_size * 2 - n_dim * window_size)
  bit_sz <- bit_sz + n_dim * window_size * log2(n_val) + n_val * n_bit

  return(list(bit_sz = bit_sz, dim_id = dim_id))
}

#' Discretize a time series using split points.
#'
#' @param motif a `matrix` with the motif.
#' @param split_pt split points for discretization.
#'
#' @return Returns the discretized time series.
#' @keywords internal
#' @noRd
#'
discretization <- function(motif, split_pt) {
  if (is.vector(motif)) {
    motif <- as.matrix(motif)
  }

  dimmotif <- dim(motif)

  for (i in 1:dimmotif[2]) {
    motif[, i] <- (motif[, i] - mean(motif[, i])) / std(motif[, i])
  }

  disc <- matrix(0, dimmotif[1], dimmotif[2])

  for (i in seq_len(length(split_pt))) {
    disc[motif < split_pt[i] & disc == 0] <- i
  }

  disc[disc == 0] <- length(split_pt) + 1

  return(disc)
}

#' Get the split points for discretization
#'
#' @param n_bit number of bits for discretization.
#'
#' @return Returns the split points.
#' @keywords internal
#' @noRd
#'
get_desc_split_pt <- function(n_bit) {
  split_pt <- stats::qnorm((1:((2^n_bit) - 1)) / (2^n_bit), 0, 1)
  return(split_pt)
}

# Global constants ---------------------------------------------------------------------------

#' Global constants
#'
#' @return Returns a `list` with the global constants
#' @keywords internal
#' @noRd
#'
vars <- function() {
  eps <- .Machine$double.eps^0.5

  return(list(eps = eps))
}

# Misc -------------------------------------------------------------------------------------------

#' Add class on front or move it to front if already exists
#'
#' @param classes result from `class()`
#' @param new_class string with the new class to add or update (put on front).
#'
#' @return Returns a vector with classes names.
#' @keywords internal
#' @noRd
#'

update_class <- function(classes, new_class) {
  classes <- classes[!(classes == new_class)]
  classes <- c(new_class, classes)

  return(classes)
}

#' Convert a TSMP object into another if possible
#'
#' The base Classes are `MatrixProfile` and `MultiMatrixProfile`, but as other functions are used,
#' classes are pushed behind, since the last output normally is the most significant. If you want,
#' for example, to plot the Matrix Profile from a `Fluss` object, you may use `as.matrixprofile()`
#' to cast it back.
#'
#' @param .mp a TSMP object.
#'
#' @describeIn as.matrixprofile Cast an object changed by another function back to `MatrixProfile`.
#' @export
#' @examples
#' \dontrun{
#'   plot(as.matrixprofile(fluss_obj))
#' }
#'

as.matrixprofile <- function(.mp) {
  if (!any(class(.mp) %in% c("MatrixProfile"))) {
    stop("Error: This object cannot be a `MatrixProfile`.")
  }

  class(.mp) <- update_class(class(.mp), "MatrixProfile")

  return(.mp)
}

#' @describeIn as.matrixprofile Cast an object changed by another function back to `MultiMatrixProfile`.
#' @export
#'

as.multimatrixprofile <- function(.mp) {
  if (!any(class(.mp) %in% c("MultiMatrixProfile"))) {
    stop("Error: This object cannot be a `MultiMatrixProfile`.")
  }

  class(.mp) <- update_class(class(.mp), "MultiMatrixProfile")

  return(.mp)
}


#' @describeIn as.matrixprofile Cast an object changed by another function back to `Fluss`.
#' @export

as.fluss <- function(.mp) {
  if (!any(class(.mp) %in% c("Fluss"))) {
    stop("Error: This object cannot be a `Fluss`.")
  }

  class(.mp) <- update_class(class(.mp), "Fluss")

  return(.mp)
}


#' @describeIn as.matrixprofile Cast an object changed by another function back to `Chain`.
#' @export

as.chain <- function(.mp) {
  if (!any(class(.mp) %in% c("Chain"))) {
    stop("Error: This object cannot be a `Chain`.")
  }

  class(.mp) <- update_class(class(.mp), "Chain")

  return(.mp)
}


#' @describeIn as.matrixprofile Cast an object changed by another function back to `Motif`.
#' @export

as.motif <- function(.mp) {
  if (!any(class(.mp) %in% c("Motif"))) {
    stop("Error: This object cannot be a `Motif`.")
  }

  class(.mp) <- update_class(class(.mp), "Motif")

  return(.mp)
}


#' @describeIn as.matrixprofile Cast an object changed by another function back to `MultiMotif`.
#' @export

as.multimotif <- function(.mp) {
  if (!any(class(.mp) %in% c("MultiMotif"))) {
    stop("Error: This object cannot be a `MultiMotif`.")
  }

  class(.mp) <- update_class(class(.mp), "MultiMotif")

  return(.mp)
}


#' @describeIn as.matrixprofile Cast an object changed by another function back to `ArcCount`.
#' @export

as.arccount <- function(.mp) {
  if (!any(class(.mp) %in% c("ArcCount"))) {
    stop("Error: This object cannot be a `ArcCount`.")
  }

  class(.mp) <- update_class(class(.mp), "ArcCount")

  return(.mp)
}

#' @describeIn as.matrixprofile Cast an object changed by another function back to `Salient`.
#' @export

as.salient <- function(.mp) {
  if (!any(class(.mp) %in% c("Salient"))) {
    stop("Error: This object cannot be a `Salient`.")
  }

  class(.mp) <- update_class(class(.mp), "Salient")

  return(.mp)
}

#' Play sound with `audio`
#'
#' @param data sound data provided by this package
#'
#' @keywords internal
#' @noRd
#'
beep <- function(data) {
  if (!(is.null(audio::audio.drivers()) || nrow(audio::audio.drivers()) == 0)) {
    tryCatch({
      audio::play(data)
    },
    error = function(cond) {
      message("Warning: Failed to play audio alert")
      message(cond)
    },
    warning = function(cond) {
      message("Warning: Something went wrong playing audio alert")
      message(cond)
    }
    )
  }
  Sys.sleep(1)
  invisible()
}
