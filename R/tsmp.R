#' Computes the Matrix Profile and Profile Index
#'
#' This is a wrap function that makes easy to use all available algorithms to compute the Matrix
#' Profile and Profile Index for multiple purposes.
#'
#' @details The Matrix Profile, has the potential to revolutionize time series data mining because
#'   of its generality, versatility, simplicity and scalability. In particular it has implications
#'   for time series motif discovery, time series joins, shapelet discovery (classification),
#'   density estimation, semantic segmentation, visualization, rule discovery, clustering etc.
#'
#'   The first algorithm invented was the [stamp()] that using [mass()] as an ultra-fast Algorithm
#'   for Similarity Search allowed to compute the Matrix Profile in reasonable time. One of its main
#'   feature was its Anytime property which using a randomized approach could return a "best-so-far"
#'   matrix that could give us the correct answer (using for example 1/10 of all iterations) almost
#'   every time.
#'
#'   The next algorithm was [stomp()] that currently is the most used. Researchers noticed that the
#'   dot products do not need to be recalculated from scratch for each subsequence. Instead, we can
#'   reuse the values calculated for the first subsequence to make a faster calculation in the next
#'   iterations. The idea is to make use of the intersections between the required products in
#'   consecutive iterations. This approach reduced the time to compute the Matrix Profile to about
#'   3% compared to [stamp()], but on the other hand, we lost the Anytime property.
#'
#'   Currently there is a new algorithm that I'll not explain further here. It is called [scrimp()],
#'   and is as fast as [stomp()], and have the Anytime property. This algorithm is implemented in
#'   this package, but still waiting for an article publication.
#'
#'   Further, there is the [mstomp()] that computes a multidimensional Matrix Profile that allows to
#'   meaningful MOTIF discovery in Multivariate Time Series. And [simple_fast()] that also handles
#'   Multivariate Time Series, but focused in Music Analysis and Exploration.
#'
#'   The [valmod()] uses a new pruning algorithm allowing a similarity search with a range of sliding
#'   window sizes.
#'
#'   The [pmp()] is a new concept that creates several profiles from a range of windows.
#'
#'   Some parameters are global across the algorithms:
#'   \describe{
#'     \item{...}{One or two time series (except for [mstomp()]). The second time series can be smaller than the first.}
#'     \item{window_size}{The sliding window.}
#'     \item{exclusion_zone}{Is used to avoid trivial matches; if a query data is provided
#'       (join similarity), this parameter is ignored.}
#'     \item{verbose}{Changes how much information is printed by this function; `0` means nothing,
#'     `1` means text, `2` adds the progress bar, `3` adds the finish sound.}
#'     \item{n_workers}{number of threads for parallel computing (except `simple_fast`, `scrimp` and `valmod`).
#'     If the value is 2 or more, the '_par' version of the algorithm will be used.}
#'   }
#'
#'   `s_size` is used only in Anytime algorithms: [stamp()] and [scrimp()].
#'   `must_dim` and `exc_dim` are used only in [mstomp()].
#'   `heap_size` is used only for [valmod()]
#'   `mode` can be any of the following: `stomp`, `stamp`, `simple`, `mstomp`, `scrimp`, `valmod`, `pmp`.
#'
#' @param \dots a `matrix` or a `vector`. If a second time series is supplied it will be a join matrix
#'   profile (except for [mstomp()]).
#' @param window_size an `int` with the size of the sliding window. Use a vector for Valmod.
#' @param exclusion_zone a `numeric`. Size of the exclusion zone, based on window size (default is
#'   `1/2`). See details.
#' @param verbose an `int`. (Default is `2`). See details.
#' @param n_workers an `int`. Number of workers for parallel. (Default is `1`).
#' @param mode the algorithm that will be used to compute the matrix profile. (Default is `stomp`).
#'   See details.
#' @param s_size a `numeric`. for anytime algorithm, represents the size (in observations) the
#'   random calculation will occur (default is `Inf`). See details.
#' @param must_dim an `int` or `vector` of which dimensions to forcibly include (default is `NULL`).
#'   See details.
#' @param exc_dim an `int` or `vector` of which dimensions to exclude (default is `NULL`). See
#'   details.
#' @param heap_size an `int`. (Default is `50`). Size of the distance profile heap buffer.
#' @param paa an `int`. (Default is `1`). Factor of PAA reduction (2 == half of size)
#' @param .keep_data a `logical`. (Default is `TRUE`). Keeps the data embedded to resultant object.
#'
#' @return Returns the matrix profile `mp` and profile index `pi`. It also returns the left and
#'   right matrix profile `lmp`, `rmp` and profile index `lpi`, `rpi` that may be used to detect
#'   Time Series Chains. [mstomp()] returns a multidimensional Matrix Profile.
#' @export
#' @references * Silva D, Yeh C, Batista G, Keogh E. Simple: Assessing Music Similarity Using
#'   Subsequences Joins. Proc 17th ISMIR Conf. 2016;23-30.
#' @references * Silva DF, Yeh C-CM, Zhu Y, Batista G, Keogh E. Fast Similarity Matrix Profile for
#'   Music Analysis and Exploration. IEEE Trans Multimed. 2018;14(8):1-1.
#' @references * Yeh CM, Kavantzas N, Keogh E. Matrix Profile VI : Meaningful Multidimensional Motif
#'   Discovery.
#' @references * Yeh CCM, Zhu Y, Ulanova L, Begum N, Ding Y, Dau HA, et al. Matrix profile I: All
#'   pairs similarity joins for time series: A unifying view that includes motifs, discords and
#'   shapelets. Proc - IEEE Int Conf Data Mining, ICDM. 2017;1317-22.
#' @references * Zhu Y, Imamura M, Nikovski D, Keogh E. Matrix Profile VII: Time Series Chains: A
#'   New Primitive for Time Series Data Mining. Knowl Inf Syst. 2018 Jun 2;1-27.
#' @references * Zhu Y, Zimmerman Z, Senobari NS, Yeh CM, Funning G. Matrix Profile II : Exploiting
#'   a Novel Algorithm and GPUs to Break the One Hundred Million Barrier for Time Series Motifs and
#'   Joins. Icdm. 2016 Jan 22;54(1):739-48.
#' @references Website: <https://sites.google.com/view/simple-fast>
#' @references Website: <https://sites.google.com/site/ismir2016simple/home>
#' @references Website: <http://www.cs.ucr.edu/~eamonn/MatrixProfile.html>
#' @family matrix profile computations
#' @examples
#' # default with [stomp()]
#' mp <- tsmp(mp_toy_data$data[1:200, 1], window_size = 30, verbose = 0)
#'
#' # Anytime STAMP
#' mp <- tsmp(mp_toy_data$data[1:200, 1], window_size = 30, mode = "stamp", s_size = 50, verbose = 0)
#'
#' # [mstomp()]
#' mp <- tsmp(mp_toy_data$data[1:200, ], window_size = 30, mode = "mstomp", verbose = 0)
#'
#' # [simple_fast()]
#' mp <- tsmp(mp_toy_data$data[1:200, ], window_size = 30, mode = "simple", verbose = 0)
#' \dontrun{
#' # parallel with [stomp_par()]
#' mp <- tsmp(mp_test_data$train$data[1:1000, 1], window_size = 30, n_workers = 2, verbose = 0)
#' }
tsmp <- function(..., window_size, exclusion_zone = getOption("tsmp.exclusion_zone", 1 / 2),
                 mode = c("stomp", "stamp", "simple", "mstomp", "scrimp", "valmod", "pmp"),
                 verbose = getOption("tsmp.verbose", 2), n_workers = 1, s_size = Inf, must_dim = NULL, exc_dim = NULL,
                 heap_size = 50, paa = 1, .keep_data = TRUE) {
  algo <- match.arg(mode)

  argv <- list(...)
  argc <- length(argv)

  if (argc == 0) {
    stop("You must supply at least one time series.")
  }

  if (argc == 1) {
    data <- argv[[1]]
    query <- NULL
  } else {
    if (argc > 2) {
      warning("Warning: Only the first two time series will be used.")
    }

    data <- argv[[1]]
    query <- argv[[2]]
  }

  paa <- round(paa)

  if (paa > 1) {
    if (is.matrix(data)) {
      data <- apply(data, 2, paa, paa)
    } else {
      data <- paa(data, paa)
    }

    if (!is.null(query)) {
      if (is.matrix(query)) {
        query <- apply(query, 2, paa, paa)
      } else {
        query <- paa(query, paa)
      }
    }

    window_size <- window_size / paa
  }

  if (n_workers > 1) {
    min_size <- length(data)

    if (min_size < 1000) {
      message("Notice: data is smaller than 1000. Single-thread mode will be used.")
      n_workers <- 1
    }
  }

  data <- as.matrix(data)
  query <- if (is.null(query)) NULL else as.matrix(query)

  result <- switch(algo,
    "stomp" = {
      if (n_workers > 1) {
        stomp_par(data, query,
          window_size = min(window_size), exclusion_zone = exclusion_zone,
          verbose = verbose, n_workers = n_workers
        )
      } else {
        stomp(data, query,
          window_size = min(window_size), exclusion_zone = exclusion_zone,
          verbose = verbose
        )
      }
    },
    "stamp" = {
      if (n_workers > 1) {
        stamp_par(data, query,
          window_size = min(window_size), exclusion_zone = exclusion_zone,
          verbose = verbose, s_size = s_size, n_workers = n_workers
        )
      } else {
        stamp(data, query,
          window_size = min(window_size), exclusion_zone = exclusion_zone,
          verbose = verbose, s_size = s_size
        )
      }
    },
    "simple" = {
      simple_fast(data, query,
        window_size = min(window_size), exclusion_zone = exclusion_zone,
        verbose = verbose
      )
    },
    "mstomp" = {
      if (argc > 1) {
        warning("Warning: Only the first time series will be used in `mstomp`.")
      }

      if (n_workers > 1) {
        mstomp_par(data,
          window_size = min(window_size), exclusion_zone = exclusion_zone,
          verbose = verbose, must_dim = must_dim, exc_dim = exc_dim, n_workers = n_workers
        )
      } else {
        mstomp(data,
          window_size = min(window_size), exclusion_zone = exclusion_zone,
          verbose = verbose, must_dim = must_dim, exc_dim = exc_dim
        )
      }
    },
    "scrimp" = {
      scrimp(data, query,
        window_size = min(window_size), exclusion_zone = exclusion_zone,
        verbose = verbose, s_size = s_size
      )
    },
    "valmod" = {
      valmod(data, query,
        window_min = min(window_size), window_max = max(window_size), heap_size = heap_size, exclusion_zone = exclusion_zone,
        verbose = verbose
      )
    },
    "pmp" = {
      pmp(data, window_sizes = window_size, n_workers = n_workers, verbose = verbose)
    },
    stop("`mode` must be ", mode)
  )

  # if (paa > 1) {
  #   result$mp <- ipaa(result$mp * sqrt(paa), paa)
  #   result$rmp <- ipaa(result$rmp * sqrt(paa), paa)
  #   result$lmp <- ipaa(result$lmp * sqrt(paa), paa)
  #   result$pi <- ipaa(result$pi, paa) * paa
  #   result$rpi <- ipaa(result$rpi, paa) * paa
  #   result$lpi <- ipaa(result$lpi, paa) * paa
  #   result$w <- result$w * paa
  #   result$paa <- paa
  #
  #   if (is.matrix(data)) {
  #     data <- apply(data, 2, ipaa, paa)
  #   } else {
  #     data <- ipaa(data, paa)
  #   }
  #
  #   if (!is.null(query)) {
  #     if (is.matrix(query)) {
  #       query <- apply(query, 2, ipaa, paa)
  #     } else {
  #       query <- ipaa(query, paa)
  #     }
  #   }
  # }

  attr(result, "origin") <- list(
    data_size = nrow(data),
    query_size = nrow(query),
    window_size = window_size,
    exclusion_zone = result$ez,
    mp_size = nrow(result$mp),
    algorithm = algo,
    class = class(result),
    version = 1.1
  )

  if (.keep_data) {
    if (!is.null(query)) {
      result$data <- list(data, query)
    } else {
      result$data <- list(data)
    }
    result$data <- lapply(result$data, as.matrix)
  }

  return(result)
}
