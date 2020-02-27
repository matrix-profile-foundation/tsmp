#' Pan-Matrix Profile
#'
#' Computes the Pan-Matrix Profile (PMP) for the given time series.
#'
#' The work closest in spirit to ours is VALMOD. The idea of VALMOD is to compute the MP for
#' the shortest length of interest, then use the information gleaned from it to guide a search
#' through longer subsequence lengths, exploiting lower bounds to prune off some calculations.
#' This idea works well for the first few of the longer subsequence lengths, but the lower bounds
#' progressively weaken, making the pruning ineffective. Thus, in the five case studies they
#' presented, the mean value of U/L was just 1.24. In contrast, consider that our termite example
#' in Fig. 15 has a U/L ratio of 240, more than two orders of magnitude larger. Thus, VALMOD is
#' perhaps best seen as finding motifs with some tolerance for a slightly (~25%) too short
#' user-specified query length, rather than a true "motif-of-all-lengths" algorithm. Also note
#' that apart from the shortest length, VALMOD only gives some information for the other lengths,
#' unlike pmp, which contains exact distances for all subsequences of all lengths.
#'
#' @param data a `matrix` or a `vector` of `numeric`.
#' @param window_sizes a `vector` of the window sizes that will be evaluated. They will be rounded to the lower integer
#' and sorted. (Default is a sequence of 20 values from 10 to half data size).
#' @param plot a `logical`. If `TRUE`, every new computation will be plotted. (Default is `FALSE`).
#' @param pmp_obj a `PMP` object that may or not contain an upper bound value, and previous computed profiles. The function will
#' add new profiles, not replace. (Default is `NULL`).
#' @param n_workers an `int`. Number of workers for parallel. (Default is `1`).
#' @param verbose an `int`. See details. (Default is `2`).
#'
#' @details
#' When just the `data` is provided, the exploration will be done using the default `window_sizes` that is a sequence
#' of 20 values between 10 and the half data size and the resulting object will have an `upper_bound` equals to `Inf`.
#' If an object is provided by the argument `pmp_obj`, this function will add more information to the resulting object,
#' never changing the values already computed.
#' `verbose` changes how much information is printed by this function; `0` means nothing, `1` means text, `2`
#' adds the progress bar, `3` adds the finish sound.
#'
#' Talk about upper bound and window sizes
#' 1. upper_window will be set to Inf on new objects
#' 1.1. upper_window will also be used for plot, and for discovery, it must not remove any existing data from the object
#' 2. window_sizes is used for plot, it must not remove any mp inside the object
#' 2.1. window_sizes tells the function what mp are stored, it may be updated with as.numeric(names(pmp))
#' 3. the functions must be capable to handle the data without need to sort by window_size, but sort may be useful later(?)
#'
#' @return Returns a `PMP` object.
#' @export
#'
#' @examples
#' # Just compute
#' pan <- pmp(mp_gait_data)
#' \dontrun{
#' # Compute the upper bound, than add new profiles
#' pan <- pmp_upper_bound(mp_gait_data)
#' pan <- pmp(mp_gait_data, pmp_obj = pan)
#' }
pmp <- function(data,
                window_sizes = seq.int(from = 10, to = length(data) / 2, length.out = 20),
                plot = FALSE,
                pmp_obj = NULL,
                n_workers = 1,
                verbose = getOption("tsmp.verbose", 2)) {

  # Parse arguments ---------------------------------
  window_sizes <- floor(window_sizes)

  checkmate::qassert(data, "N+")
  checkmate::qassert(window_sizes, "X+")
  checkmate::qassert(plot, "B")
  checkmate::qassert(pmp_obj, c("0", "l"))
  checkmate::qassert(n_workers, paste0("X1[1,", parallel::detectCores(), "]"))
  checkmate::qassert(verbose, "X1")

  # checks if the given object is actualy a skimp object
  if (!is.null(pmp_obj)) {
    if (class(pmp_obj) != "PMP") {
      stop("`pmp_obj` must be of class `PMP`")
    }
  }

  ## Prepare things ----

  data_size <- length(data)

  # if an object is given, remove the windows that already have been computed and are below the upper_window
  if (!is.null(pmp_obj)) {
    # remove already computed
    window_sizes <- window_sizes[!(window_sizes %in% pmp_obj$w)]

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

  # print(str(list(min_window = min_window, max_window = max_window, window_sizes = window_sizes)))

  # if we'll plot while computing, prepare the canvas
  if (plot == TRUE) {
    # if an object is given, plot it using existing windows
    if (!is.null(pmp_obj)) {
      graphics::plot(pmp_obj) # skimp_plot_set_canvas(pmp_obj = pmp_obj)
      Sys.sleep(1) # needed for plot update
    } else {
      # create a blank canvas with the proper size
      skimp_plot_set_canvas(
        ymin = min_window,
        ymax = max_window + floor((max_window - min_window) / 24), # arbitrary
        xmin = 1,
        xmax = data_size
      )
      Sys.sleep(1) # needed for plot update
    }
  }

  # anytime must return the result always
  on.exit(
    {
      if (is.null(pmp_obj)) {
        return(NULL)
      } else {

        # final arrangements to existing object

        expected_profiles <- as.numeric(names(pmp_obj$pmp))
        expected_indexes <- as.numeric(names(pmp_obj$pmpi))

        # check if the windows vector contains any uncomputed matrix profile
        if (any(!(pmp_obj$w %in% expected_profiles))) {
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

        pmp_obj$ez <- getOption("tsmp.exclusion_zone", 1 / 2)

        return(pmp_obj)
      }
    },
    TRUE
  )

  # if not given, create a new object to start with
  if (is.null(pmp_obj)) {
    pmp_obj <- list(pmp = list(), pmpi = list())
    class(pmp_obj) <- "PMP"
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
    result <- mpx(data = data, window_size = w, idx = TRUE, dist = "euclidean", n_workers = n_workers)

    if (verbose > 0) {
      message(
        "step: ", i, "/", length(split_idx), " binary idx: ", idx, " window: ", w
      )
    }

    # if pmp_obj is a new empty object, accessing windows will return NULL, so it's fine
    pmp_obj$w <- c(pmp_obj$w, w)
    # using character to create a tuple list. Numbers would create NULL's
    pmp_obj$pmp[[as.character(w)]] <- result$mp
    pmp_obj$pmpi[[as.character(w)]] <- result$pi

    if (plot == TRUE) {
      # add a layer to the plot. `pmp_obj$w` is currently used to know the heigth of
      # the new layer. May be room to improve.
      skimp_plot_add_layer(result$mp, w, pmp_obj$w)
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

#' Pan Matrix Profile upper bound
#'
#' Finds the upper bound for Pan Matrix Profile calculation.
#'
#' @param data a `matrix` or a `vector` of `numeric`.
#' @param threshold a `numeric`. Correlation threshold. See details.  (Default is `0.95`).
#' @param refine_stepsize a `numeric`. Step size for the last upper bound search. See details.  (Default is `0.25`).
#' @param return_pmp a `logical`. If `TRUE`, returns the computed data as a `PMP` object, if `FALSE`,
#' returns just the upper bound value. (Default is `TRUE`).
#' @param n_workers an `int`. Number of workers for parallel. (Default is `1`).
#' @param verbose verbose an `int`. See details. (Default is `2`).
#'
#' @details
#' The Pan Matrix Profile may not give any further information beyond a certain window size. This function starts
#' computing the matrix profile for the window size of 8 and doubles it until the minimum correlation value found is
#' less than the `threshold`. After that, it begins to refine the upper bound using the `refine_stepsize` values, until
#' the `threshold` value is hit.
#'
#' `verbose` changes how much information is printed by this function; `0` means nothing, `1` means text, `2`
#' adds the progress bar, `3` adds the finish sound.
#'
#' @return Returns a `PMP` object with computed data, or just the upper bound value if `return_pmp` is set to `FALSE`.
#' @export
#'
#' @references * Yet to be announced
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#'
#' @examples
#' # return the object
#' pan_matrix <- pmp_upper_bound(mp_gait_data)
#'
#' # just the upper bound
#' pan_ub <- pmp_upper_bound(mp_gait_data, return_pmp = FALSE)
pmp_upper_bound <- function(data,
                            threshold = getOption("tsmp.pmp_ub", 0.95),
                            refine_stepsize = getOption("tsmp.pmp_refine", 0.25),
                            return_pmp = TRUE, n_workers = 1,
                            verbose = getOption("tsmp.verbose", 2)) {
  correlation_max <- Inf

  if (return_pmp) {
    pmp <- list()
    pmpi <- list()
    do_idxs <- TRUE
  } else {
    do_idxs <- FALSE
  }
  correlation_max <- Inf
  window_size <- 8
  windows <- NULL
  max_window <- floor(length(data) / 2)

  # Register the anytime exit point
  on.exit(
    {
      if (return_pmp) {
        pmp_obj <- list(upper_window = max(windows), pmp = pmp, pmpi = pmpi, w = windows)
        class(pmp_obj) <- "PMP"
        return(pmp_obj)
      } else {
        return(window_size)
      }
    },
    TRUE
  )

  # first perform a wide search by increasing window by 2 in each iteration
  while (window_size <= max_window) {
    # message("window: ", window_size)

    result <- mpx(data = data, window_size = window_size, idx = do_idxs, dist = "pearson", n_workers = n_workers)
    if (is.null(result) || result$partial) {
      warning("The computation was terminated prematurely. The results are partial.")
      return()
    }

    correlation_max <- max(result$mp[!is.infinite(result$mp)], na.rm = TRUE)

    if (correlation_max < threshold) {
      # message("break at ", window_size)
      break
    }

    if (return_pmp) {
      pmp[[as.character(window_size)]] <- corr_ed(result$mp, window_size)
      pmpi[[as.character(window_size)]] <- result$pi
    }

    windows <- c(windows, window_size)

    window_size <- window_size * 2
  }

  if (window_size <= max_window) {
    # refine the upper bound by increase by + X% increments
    test_windows <- 2 * round(((seq(refine_stepsize, 1 - 1e-5, refine_stepsize) + 1) * window_size / 2) / 2)

    for (window_size in test_windows) {
      # message("refine window: ", window_size)

      result <- mpx(data = data, window_size = window_size, idx = do_idxs, dist = "pearson", n_workers = n_workers)
      if (is.null(result) || result$partial) {
        warning("The computation was terminated prematurely. The results are partial.")
        return()
      }

      windows <- c(windows, window_size)

      correlation_max <- max(result$mp[!is.infinite(result$mp)], na.rm = TRUE)

      if (return_pmp) {
        pmp[[as.character(window_size)]] <- corr_ed(result$mp, window_size)
        pmpi[[as.character(window_size)]] <- result$pi
      }

      if (correlation_max < threshold) {
        # message("break refine at ", window_size)
        break
      }
    }
  }
}
