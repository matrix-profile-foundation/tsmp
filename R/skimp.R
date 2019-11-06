#' Title
#'
#' @param data dataframe
#' @param window_sizes window sizes
#'
#'   The work closest in spirit to ours is VALMOD. The idea of VALMOD is to compute the MP for
#'   the shortest length of interest, then use the information gleaned from it to guide a search
#'   through longer subsequence lengths, exploiting lower bounds to prune off some calculations.
#'   This idea works well for the first few of the longer subsequence lengths, but the lower bounds
#'   progressively weaken, making the pruning ineffective. Thus, in the five case studies they
#'   presented, the mean value of U/L was just 1.24. In contrast, consider that our termite example
#'   in Fig. 15 has a U/L ratio of 240, more than two orders of magnitude larger. Thus, VALMOD is
#'   perhaps best seen as finding motifs with some tolerance for a slightly (~25%) too short
#'   user-specified query length, rather than a true "motif-of-all-lengths" algorithm. Also note
#'   that apart from the shortest length, VALMOD only gives some information for the other lengths,
#'   unlike pmp, which contains exact distances for all subsequences of all lengths.
#'
#' @param plot plot each interaction
#' @param pmp_res
#'
#' @return Returns the Pan Matrix Profile
#' @export
#'
#' @examples
skimp <- function(data,
                  window_sizes = seq.int(from = 10, to = length(data) / 2, length.out = 20),
                  plot = FALSE,
                  pmp_obj = NULL) {
  window_sizes <- floor(window_sizes)
  data_size <- length(data)

  if (is.list(pmp_obj)) {
    window_sizes <- window_sizes[!(window_sizes %in% pmp_obj$windows)]
  } else {
    pmp_obj <- list(pmp = list(), pmpi = list(), windows = 0)
  }

  window_sizes <- sort(window_sizes)
  min_window <- min(window_sizes)
  max_window <- max(window_sizes)

  # Determine the order in which we will explore the subsequence lengths used to run the matrix profile
  split_idx <- binary_split(length(window_sizes))
  #split_idx <- seq_along(window_sizes)


  ## BEGIN PROCESSING

  # Runs matrix profile on each subsequence length in INPUTS[['window_sizes']](split_idx), stores the
  # results in pmp, and generates a frame for the animation

  plot(c(1, data_size), c(min_window, max_window + diff(window_sizes[1:2])),
    main = "Pan Matrix Profile",
    type = "n", xlab = "", ylab = "window"
  )

  for (i in 1:length(split_idx)) {

    # Get the current subsequence length to run the Matrix Profile on
    w <- window_sizes[split_idx[i]]

    if (w == 0) {
      next
    }

    # Run Matrix Profile
    result <- mpx(data, w, idx = TRUE, dist = "euclidean")

    message(
      "w: ", w, " i: ", i, "/", length(split_idx), " idx: ", split_idx[i]
    )

    pmp_obj$pmp[[as.character(w)]] <- result$mp
    pmp_obj$pmpi[[as.character(w)]] <- result$pi

    test <- c(result$mp, rep(0, w - min_window))
    test[test > 1] <- 1
    test[test < 0] <- 0

    heigth <- window_sizes[split_idx[1:i]]
    heigth <- sort(heigth[heigth > w])[1]
    if (is.na(heigth)) {
      heigth <- max_window + diff(window_sizes[1:2])
    }

    message("from: ", w, " to: ", heigth)

    graphics::rasterImage(matrix(test, nrow = 1),
      xleft = 1, xright = length(test),
      ybottom = w - 0.2, ytop = heigth,
      interpolate = FALSE
    )
    dev.flush()
  }

  sorted <- sort(as.numeric(names(pmp_obj$pmp)), index.return = TRUE)
  sort_idxs <- sorted$ix
  window_sizes <- sorted$x

  pmp_obj$pmp <- pmp_obj$pmp[sort_idxs]
  pmp_obj$pmpi <- pmp_obj$pmpi[sort_idxs]
  pmp_obj$windows <- window_sizes

  return(pmp_obj)
  #   A dict with the following:
  #   {
  #     'pmp': the pan matrix profile as a 2D array,
  #     'pmpi': the pmp indices,
  #     'data': {
  #       'ts': time series used,
  #     },
  #     'windows': the windows used to compute the pmp,
  #     'sample_pct': the sample percent used,
  #     'metric':The distance metric computed for the pmp,
  #     'algorithm': the algorithm used,
  #     'class': PMP
  #   }
} # function skimp

#' Title
#'
#' @param pmp
#' @param inv
#'
#' @return
#' @export
#'
#' @examples
plot_skimp <- function(pmp, cr = 0, inv = FALSE) {
  test <- list_to_matrix(pmp$pmp)
  xmax <- ncol(test)
  ymin <- min(pmp$windows)
  ymax <- max(pmp$windows)
  sizes <- diff(pmp$windows)
  sizes <- c(sizes, tail(sizes, 1))
  data_size <- xmax + ymin - 1

  # depth <- 256
  # test <- ceiling(test * depth) / depth
  #
  # m <- max(test[!is.infinite(test)])
  # n <- min(test)
  # test <- (test - n) / (m - n)

  test[test > 1] <- 1
  test[test < 0] <- 0

  if (inv) {
    test <- 1 - test
  }

  plot(c(1, data_size), c(ymin, ymax + diff(pmp$windows[1:2])),
    main = "Pan Matrix Profile",
    type = "n", xlab = "", ylab = "window"
  )

  for (i in seq_along(pmp$windows)) {
    graphics::rasterImage(test[i, , drop = FALSE],
      xleft = 1, xright = xmax,
      ybottom = pmp$windows[i], ytop = pmp$windows[i] + sizes[i],
      interpolate = FALSE
    )
  }

  dev.flush()
}
# plot_skimp <- function(pmp, inv = FALSE) {
#   test <- pmp$pmp[rev(seq_along(pmp$windows)), ]
#
#   depth <- 256
#   test <- ceiling(test * depth) / depth
#
#   # m <- max(test[!is.infinite(test)])
#   # n <- min(test)
#   # test <- (test - n) / (m - n)
#
#   test[test > 1] <- 1
#   test[test < 0] <- 0
#
#   if (inv) {
#     test <- 1 - test
#   }
#
#   plot.new()
#   grid::grid.raster(test, height = 0.5, width = 0.8, interpolate = FALSE)
# }

#' Pan Matrix Profile maximum subsequence
#'
#' Finds the upper bound
#'
#' @param data
#' @param threshold
#' @param refine_stepsize
#' @param return_pmp
#'
#' @return
#' @export
#'
#' @examples
#'
maximum_subsequence <- function(data,
                                threshold = getOption("tsmp.pmp_ub", 0.95),
                                refine_stepsize = getOption("tsmp.pmp_refine", 0.5),
                                return_pmp = TRUE) {
  correlation_max <- Inf

  if (return_pmp) {
    pmp <- list()
    pmpi <- list()
    do_idxs <- TRUE
  }
  correlation_max <- Inf
  window_size <- 8
  windows <- NULL
  max_window <- floor(length(data) / 2)

  # first perform a wide search by increasing window by 2 in each iteration
  while (window_size <= max_window) {
    message("window: ", window_size)

    windows <- c(windows, window_size)

    result <- mpx(data, window_size, idx = do_idxs, dist = "pearson")
    correlation_max <- max(result$mp[!is.infinite(result$mp)], na.rm = TRUE)

    if (return_pmp) {
      pmp[[as.character(window_size)]] <- corr_ed(result$mp, window_size)
      pmpi[[as.character(window_size)]] <- result$pi
    }

    if (correlation_max < threshold) {
      break
    }

    window_size <- window_size * 2
  }

  # refine the upper bound by increase by + X% increments
  test_windows <- 2 * round(((seq(refine_stepsize, 1 - refine_stepsize, refine_stepsize) + 1) * window_size / 2) / 2)

  for (window_size in test_windows) {
    message("window: ", window_size)

    windows <- c(windows, window_size)

    result <- mpx(data, window_size, idx = do_idxs, dist = "pearson")
    correlation_max <- max(result$mp[!is.infinite(result$mp)], na.rm = TRUE)

    if (return_pmp) {
      pmp[[as.character(window_size)]] <- corr_ed(result$mp, window_size)
      pmpi[[as.character(window_size)]] <- result$pi
    }

    if (correlation_max < threshold) {
      break
    }
  }

  if (return_pmp) {
    return(list(upper_window = window_size, pmp = pmp, pmpi = pmpi, windows = windows))
  } else {
    return(window_size)
  }
}
