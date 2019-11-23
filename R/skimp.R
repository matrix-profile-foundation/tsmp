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
#'
#' 1. upper_window will be set to Inf on new objects
#' 1.1. upper_window will also be used for plot, and for discovery, it must not remove any existing data from the object
#' 2. window_sizes is used for plot, it must not remove any mp inside the object
#' 2.1. window_sizes tells the function what mp are stored, it may be updated with as.numeric(names(pmp))
#' 3. the functions must be capable to handle the data without need to sort by window_size, but sort may be useful later(?)

skimp <- function(data,
                  window_sizes = seq.int(from = 10, to = length(data) / 2, length.out = 20),
                  plot = FALSE,
                  sorted = FALSE,
                  pmp_obj = NULL) {

  ## Prepare things ----

  # checks if the given object is actualy a skimp object
  if (!is.null(pmp_obj)) {
    if (class(pmp_obj) != "skimp") {
      stop("`pmp_obj` must be of class skimp")
    }
  }

  # TODO: sanitize `data`
  data_size <- length(data)

  # windows must be integers
  window_sizes <- floor(window_sizes)

  # if an object is given, remove the windows that already have been computed and are below the upper_window
  if (!is.null(pmp_obj)) {
    # remove already computed
    window_sizes <- window_sizes[!(window_sizes %in% pmp_obj$windows)]

    if (!is.null(pmp_obj$upper_window)) {
      # remove those above the upper_window
      window_sizes <- window_sizes[window_sizes < pmp_obj$upper_window]
    } else {
      # else, keep them and set the upper_window to Inf, since it is NULL
      pmp_obj$upper_window <- Inf
    }
  }

  # get the extreme values and sort the window_sizes for later binary_split correctly
  min_window <- min(window_sizes)
  max_window <- max(window_sizes)
  window_sizes <- sort(window_sizes)

  print(str(list(min_window = min_window, max_window = max_window, window_sizes = window_sizes)))

  # if we'll plot while computing, prepare the canvas
  if (plot == TRUE) {
    # if an object is given, plot it using existing windows
    if (!is.null(pmp_obj)) {
      plot(pmp_obj, pearson = T) # skimp_plot_set_canvas(pmp_obj = pmp_obj)
      Sys.sleep(1) # needed for plot update
    } else {
      # create a blank canvas with the proper size
      skimp_plot_set_canvas(ymin = min_window,
                            ymax = max_window + floor((max_window - min_window) / 24), # arbitrary
                            xmin = 1,
                            xmax = data_size)
      Sys.sleep(1) # needed for plot update
    }
  }

  # anytime must return the result always
  on.exit({

    if (is.null(pmp_obj)) {
      message("on.exit NULL")

      return(NULL)
    } else {

      # final arrangements to existing object

      expected_profiles <- as.numeric(names(pmp_obj$pmp))
      expected_indexes <- as.numeric(names(pmp_obj$pmpi))

      # check if the windows vector contains any uncomputed matrix profile
      if (any(!(pmp_obj$windows %in% expected_profiles))) {
        warning("`windows` contains values not computed in `pmp`")
      }

      # check if there is any matrix profile without a profile index
      if (any(!(expected_profiles %in% expected_indexes))) {
        warning("`pmp` contains values without respective `pmpi")
      }

      # check if there is any profile index without a matrix profile
      if (any(!(expected_indexes %in% expected_profiles))) {
        warning("`pmpi` contains values without respective `pmp")
      }


      # if (sorted == TRUE) {
      #   sorted <- sort(as.numeric(names(pmp_obj$pmp)), index.return = TRUE)
      #   window_sizes <- sorted$x
      #
      #   pmp_obj$pmp <- pmp_obj$pmp[sorted$ix]
      #   pmp_obj$pmpi <- pmp_obj$pmpi[sorted$ix]
      # } else {
      #   window_sizes <- as.numeric(names(pmp_obj$pmp))
      # }

      message("on.exit OBJ")

      # pmp_obj$windows <- window_sizes
      #
      class(pmp_obj) <- "skimp"

      return(pmp_obj)
    }
  }, TRUE)

  # if not given, create a new object to start with
  if (is.null(pmp_obj)) {
    pmp_obj <- list(pmp = list(), pmpi = list())
  }

  # Determine the order in which we will explore the window_sizes
  split_idx <- binary_split(length(window_sizes))

  ## Begin the main Loop ----

  for (i in seq_along(split_idx)) {

    # i holds the sequence from 1 to ...length(split_idx)
    # idx holds the actual binary split index
    idx <- split_idx[i]

    # w holds the current window being explored
    w <- window_sizes[idx]

    if (is.na(w) || w == 0) {
      warning("Invalid window size ", w)
      next
    }

    # Run Matrix Profile
    result <- mpx(data, w, idx = TRUE, dist = "euclidean")

    message(
      "step: ", i, "/", length(split_idx), " binary idx: ", idx, " window: ", w
    )

    # if pmp_obj is a new empty object, accessing windows will return NULL, so it's fine
    pmp_obj$windows <- c(pmp_obj$windows, w)
    # using character to create a tuple list. Numbers would create NULL's
    pmp_obj$pmp[[as.character(w)]] <- result$mp
    pmp_obj$pmpi[[as.character(w)]] <- result$pi

    if (plot == TRUE) {
      # add a layer to the plot. `pmp_obj$windows` is currently used to know the heigth of
      # the new layer. May be room to improve.
      skimp_plot_add_layer(result$mp, w, pmp_obj$windows)
      Sys.sleep(1) # needed for plot update
    }
  }


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



skimp_plot_set_canvas <- function(..., pmp_obj = NULL) {
  if (!is.null(pmp_obj)) {
    xmin <- 1
    ymin <- min(pmp_obj$windows)
    ymax <- max(pmp_obj$windows) + floor((max(pmp_obj$windows) - ymin) / 24) # arbitrary
    mp_min <- length(pmp_obj$pmp[[as.character(ymin)]])
    xmax <- mp_min + ymin - 1
  } else {
    pars <- list(...)
    xmin <- pars$xmin
    xmax <- pars$xmax
    ymin <- pars$ymin
    ymax <- pars$ymax

    plot(c(xmin, xmax), c(ymin, ymax),
         main = "Pan Matrix Profile", xlab = "", ylab = "window",
         type = "n", xaxt = "n", yaxt = "n", xlim = c(xmin, xmax)
    )

    axis(side = 1, at = floor(c(xmin, seq(xmin, xmax, length.out = 10), xmax))) # X
    axis(side = 2, at = floor(c(ymin, seq(ymin, ymax, length.out = 10), ymax))) # Y
  }

  invisible()
}

skimp_plot_add_layer <- function(layer, window, window_set = NULL) {
  coords <- par("usr")
  xmin <- 1
  ymin <- window
  data_size <- length(layer) + window - 1 # theoretical data size

  # assert
  if (data_size != (coords[[2]] + coords[[1]] - xmin)) {
    stop("data_size calc is wrong")
  }

  if (is.null(window_set)) {
    # if this is the first layer
    w_min <- window
  } else {
    # else, find the position to plot
    w_min <- min(window_set)
    ymax <- coords[[4]] + coords[[3]] - w_min
    # assert
    if (ymax != (max(window_set) + floor((max(window_set) - min(window_set)) / 24))) {
      print(list(
        ymax = ymax,
        max_window = max(window_set),
        min_window = min(window_set),
        result = (max(window_set) + floor((max(window_set) - min(window_set)) / 24))
      ))
      stop("ymax calc is wrong")
    }

    upper_windows <- window_set[window_set > window]

    if (length(upper_windows) > 0) {
      ytop <- min(upper_windows)
    } else {
      ytop <- ymax
    }
  }

  layer <- c(layer, rep(0, window - 1))
#  layer <- ed_corr(layer, window)
  xmax <- length(layer)

  layer <- normalize(layer, 0, 1)

  # layer[layer > 1] <- 1
  # layer[layer < 0] <- 0

  # print(list(
  #   coords = coords, xmin = xmin, xmax = xmax,
  #   ymin = ymin, ymax = ymax, ytop = ytop, w_min = w_min
  # ))

  graphics::rasterImage(matrix(layer, nrow = 1),
    xleft = xmin, xright = xmax,
    ybottom = ymin, ytop = ytop,
    interpolate = FALSE
  )

  invisible()
}



#' Title
#'
#' @param pmp
#' @param inv
#'
#' @return
#' @export
#'
#' @examples
plot.skimp <- function(pmp, pearson = FALSE) {
  sorted <- sort(pmp$windows, index.return = TRUE)

  min_window <- min(sorted$x)
  max_window <- max(sorted$x)
  data_size <- length(pmp$pmp[[1]]) + min_window - 1
  xmin <- 1
  xmax <- data_size - min_window + 1
  ymin <- min_window

  sizes <- diff(sorted$x)
  sizes <- c(sizes, sizes[1])

  if (length(pmp$windows) > 1) {
    ymax <- max_window + floor((max_window - min_window) / 24) # arbitrary
  } else {
    ymax <- max_window + floor((max_window - min_window) / 24) # arbitrary
  }

  plot(c(xmin, data_size), c(ymin, ymax),
    main = "Pan Matrix Profile", xlab = "", ylab = "window",
    type = "n", xaxt = "n", yaxt = "n", xlim = c(xmin, data_size)
  )

  axis(side = 1, at = floor(c(xmin, seq(xmin, data_size, length.out = 10), data_size))) # X
  axis(side = 2, at = floor(c(ymin, seq(ymin, ymax, length.out = 10), ymax))) # Y

  for (i in seq_along(sorted$ix)) {
    test <- pmp$pmp[[sorted$ix[i]]]
    w <- pmp$windows[sorted$ix[i]]

    test <- c(test, rep(0, w - min_window))

    # if (pearson == TRUE) {
    #   test <- ed_corr(test, w)
    # }

    test <- normalize(test, 0, 1)

    ytop <- w + sizes[sorted$ix[i]] * 1.5

    if (ytop > ymax) {
      ytop <- ymax
    }

    graphics::rasterImage(matrix(test, nrow = 1),
      xleft = 1, xright = xmax,
      ybottom = w, ytop = ytop,
      interpolate = FALSE
    )
  }
}

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

    result <- mpx(data, window_size, idx = do_idxs, dist = "pearson")
    correlation_max <- max(result$mp[!is.infinite(result$mp)], na.rm = TRUE)

    if (correlation_max < threshold) {
      break
    }

    if (return_pmp) {
      pmp[[as.character(window_size)]] <- corr_ed(result$mp, window_size)
      pmpi[[as.character(window_size)]] <- result$pi
    }

    windows <- c(windows, window_size)

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
    pmp_obj <- list(upper_window = window_size, pmp = pmp, pmpi = pmpi, windows = windows)
    class(pmp_obj) <- "skimp"
    return(pmp_obj)
  } else {
    return(window_size)
  }
}
