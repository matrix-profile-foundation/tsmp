#' Calculates the distance profile using MASS algorithms
#'
#' Mueen's Algorithm for Similarity Search is The Fastest Similarity Search Algorithm for Time
#' Series Subsequences under Euclidean Distance and Correlation Coefficient.
#'
#' @details
#' This function has several ways to work:
#'
#' Case 1: You have a small sized query and the data. In this case you only have to provide the first two
#' parameters `data` and `query`. Internally the `window_size` will be get from the query length.
#'
#' Case 2: You have one or two data vectors and want to compute the join or self-similarity. In this case
#' you need to use the recursive solution. The parameters are `data`, `query`, `window_size` and `index`.
#' The first iteration don't need the `index` unless you are starting somewhere else. The `query` will be
#' the source of a `query_window`, starting on `index`, with length of `window_size`.
#'
#' The `method` defines which MASS will be used. Current supported values are: `v2`, `v3`, `weighted`.
#'
#' @param data a `matrix` or a `vector`.
#' @param query a `matrix` or a `vector`. See details.
#' @param params Precomputed values from the first iteration. If not supplied, these values will be computed.
#' @param window_size an `int` or `NULL`. Sliding window size. See details.
#' @param index an `int`. Index of query window. See details.
#' @param method method that will be used to calculate the distance profile. See details.
#' @param type comparison type used by MASS: `"normalized"`, `"non_normalized"`, `"absolute"`, or `"weighted"`.
#' @param weights,weight a `vector` of `numeric` or `NULL` with the same length of the `window_size`. This is
#' a MASS extension to weight the query.
#' @param n_workers an `int`. Number of threads to use when computing with MASS. (Default is `1`).
#' @param \dots Placeholder for legacy arguments.
#'
#' @return Returns the `distance_profile` for the given query and the `last_product` for STOMP
#'   algorithm and the parameters for recursive call. See details.
#'
#' @export
#'
#' @references * Abdullah Mueen, Yan Zhu, Michael Yeh, Kaveh Kamgar, Krishnamurthy Viswanathan,
#'   Chetan Kumar Gupta and Eamonn Keogh (2015), The Fastest Similarity Search Algorithm for Time
#'   Series Subsequences under Euclidean Distance
#' @references Website: <https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html>
#'
#' @examples
#'
#' w <- mp_toy_data$sub_len
#' ref_data <- mp_toy_data$data[, 1]
#' # minimum example, data and query
#' nn <- dist_profile(ref_data, ref_data[1:w])
#' distance_profile <- sqrt(nn$distance_profile)
#'
#' # data and indexed query
#' nn <- dist_profile(ref_data, ref_data, window_size = w, index = 10)
#' distance_profile <- sqrt(nn$distance_profile)
#'
#' # recursive
#' nn <- NULL
#'
#' for (i in seq_len(10)) {
#'   nn <- dist_profile(ref_data, ref_data, nn, window_size = w, index = i)
#' }
#'
#' # weighted
#' weight <- c(rep(1, w / 3), rep(0.5, w / 3), rep(0.8, w / 3)) # just an example
#'
#' nn <- dist_profile(ref_data, ref_data,
#'   window_size = w, index = 1, method = "weighted",
#'   weight = weight
#' )
#' distance_profile <- sqrt(nn$distance_profile)
dist_profile <- function(data, query = data, params = NULL, window_size = NULL, index = 1,
                         method = c("v3", "v2", "weighted"),
                         type = c("normalized", "non_normalized", "absolute", "weighted"),
                         weights = NULL, weight = weights, n_workers = 1L, ...) {

  # Handle window_size passed as 2nd positional argument (e.g. dist_profile(data, 50))
  if (is.numeric(query) && length(query) == 1L && is.null(window_size)) {
    window_size <- as.integer(query)
    query <- data
  }

  # Check time series classes -----------------------
  data <- convert_data(data)
  query <- convert_data(query)

  if (is.null(window_size)) {
    if (is.list(params)) {
      if (!is.null(params$window_size)) {
        window_size <- params$window_size
      } else if (!is.null(params$par$window_size)) {
        window_size <- params$par$window_size
      } else if (!is.null(params$params$window_size)) {
        window_size <- params$params$window_size
      }
    }
    if (is.null(window_size)) {
      window_size <- length(query)
    }
  }

  # Parse arguments ---------------------------------
  "!!!DEBUG Parsing data"
  checkmate::qassert(data, "n+")
  window_size <- checkmate::qassert(window_size, "x+[4,)")
  checkmate::qassert(query, "n+")
  checkmate::qassert(index, glue::glue("X1[1,{length(data) - window_size + 1}]"))
  n_workers <- as.integer(checkmate::qassert(n_workers, "X+"))
  if (!missing(weight) && is.null(weights)) {
    weights <- weight
  }
  "!!!DEBUG Parsing type"
  method_in <- match.arg(method)
  if (method_in == "weighted") {
    type <- "weighted"
    method <- "v3"
  } else {
    method <- method_in
    type <- match.arg(type)
  }
  if (type == "weighted") {
    if (is.null(weights)) {
      stop("The `weights` argument must be provided.", call. = FALSE)
    }
    if (length(weights) != window_size) {
      stop("The `weights` must be the same size as the `window_size`.", call. = FALSE)
    }
    checkmate::qassert(weights, "N+")
  }

  # Register anytime exit point ----------------------

  result <- NULL

  "!DEBUG Register anytime exit point"
  on.exit(
    if (is.null(params)) {
      return(invisible(NULL))
    } else {
      return(c(result, list(par = params, params = params)))
    },
    TRUE
  )

  if (is.null(query)) {
    "!!!DEBUG query is null"
    query <- data
  }

  if (length(query) > window_size) {
    query_sub <- query[index:(window_size + index - 1)]
  } else {
    query_sub <- query
  }

  # Computation ------------------------------------
  "!DEBUG Computation"
  if (anyNA(query_sub)) {
    ## ---- Query with Gap ----
    q1 <- NULL
    q2 <- NULL

    min_idx <- min(which(is.na(query_sub))) - 1
    max_idx <- max(which(is.na(query_sub))) + 1

    if (min_idx >= 4L) {
      q1 <- query_sub[1:min_idx]
    }

    if (max_idx <= (window_size - 4)) {
      q2 <- query_sub[max_idx:window_size]
    }

    if (anyNA(q1) || anyNA(q2)) {
      stop("Querying with gap only supports one gap.")
    }

    q1_len <- length(q1)
    q2_len <- length(q2)

    if (q1_len != q2_len) {
      warning("Warning: Result may be inconsistent if the size of the queries are different.")
    }

    pre1 <- matrixprofiler::mass_pre(data, window_size = q1_len, query = q1)
    pre2 <- matrixprofiler::mass_pre(data, window_size = q2_len, query = q2)

    result1 <- matrixprofiler::mass(pre1, data, q1, 1, version = "v3", n_workers = n_workers)
    result2 <- matrixprofiler::mass(pre2, data, q2, 1, version = "v3", n_workers = n_workers)
    result1 <- sqrt(result1$distance_profile)
    result2 <- sqrt(result2$distance_profile)

    pad <- rep(Inf, max_idx - 1)
    result2 <- c(pad, result2)
    result1 <- c(result1, rep(Inf, length(result2) - length(result1)))

    result <- list(distance_profile = (result1 + result2)^2)
    params <- list(pre1 = pre1, pre2 = pre2, window_size = window_size)
    return(c(result, list(par = params, params = params)))
  }

  if (is.null(params)) {
    "!!!DEBUG params is null"
    ## First iteration with MASS ----
    tryCatch(
      {
        pre <- matrixprofiler::mass_pre(data, window_size, query, type, weights)
      },
      error = print
    )
    "!!!DEBUG params is a list with par and window_size"
    params <- pre
    params$window_size <- window_size
  } else {
    "!!!DEBUG params is not null"
    ## Following iterations with MASS ----
    if (!is.null(params$par)) {
      params <- params$par
    }
  }

  ## Do the search --------
  tryCatch(
    {
      result <- matrixprofiler::mass(params, data, query, index, version = method, n_workers = n_workers)
    },
    error = print
  )

  checkmate::qassert(params, "L+")
  checkmate::qassert(result, "L+")
}

mass_v2 <- function(query, window_size, data_fft, data_size, data_mean, data_sd, query_mean, query_sd) {
  pre_obj <- structure(
    list(
      data_size = as.integer(data_size),
      window_size = as.integer(window_size),
      data_fft = data_fft,
      data_mean = as.numeric(data_mean),
      data_sd = as.numeric(data_sd),
      query_mean = as.numeric(query_mean)[1],
      query_sd = as.numeric(query_sd)[1]
    ),
    type = "normalized"
  )
  dummy_data <- rep(0, data_size)
  matrixprofiler::mass(pre_obj, data = dummy_data, query = query, index = 1, version = "v2")
}

mass_v3 <- function(query, data, window_size, data_size, data_mean, data_sd, query_mean, query_sd) {
  pre_obj <- structure(
    list(
      data_size = as.integer(data_size),
      window_size = as.integer(window_size),
      data_mean = as.numeric(data_mean),
      data_sd = as.numeric(data_sd),
      query_mean = as.numeric(query_mean)[1],
      query_sd = as.numeric(query_sd)[1]
    ),
    type = "normalized"
  )
  matrixprofiler::mass(pre_obj, data = data, query = query, index = 1, version = "v3")
}

mass_pre_w <- function(data, query = NULL, window_size, weight) {
  data_size <- length(data)
  data_avgsd <- list(
    avg = matrixprofiler::mov_mean(data, window_size),
    sd = matrixprofiler::mov_std(data, window_size)
  )
  data_mean <- data_avgsd$avg
  data_sd <- data_avgsd$sd
  data_fft <- stats::fft(data)

  sumw <- sum(weight)
  rweight <- rev(weight)
  rweight[(window_size + 1):data_size] <- 0
  w_fft <- stats::fft(rweight)
  data_fft_w_fft <- data_fft * w_fft
  data_w <- stats::fft(data_fft_w_fft, inverse = TRUE) / length(data_fft_w_fft)
  data2_fft <- stats::fft(data^2)
  data2_fft_w_fft <- data2_fft * w_fft
  data2_w <- stats::fft(data2_fft_w_fft, inverse = TRUE) / length(data2_fft_w_fft)
  sumxw2 <- data2_w[window_size:data_size]
  sumxw <- data_w[window_size:data_size]

  data_pre <- ((sumxw2 - 2 * sumxw * data_mean + sumw * data_mean^2) / data_sd^2)

  return(list(
    window_size = window_size, data_fft = data_fft, data_size = data_size,
    data_mean = data_mean, data_sd = data_sd, data_pre = data_pre, weight = weight
  ))
}

mass_weighted <- function(query_window, window_size, data_fft, data_size, data_mean, data_sd, data_pre, weight, ...) {
  query_window <- (query_window - mean(query_window)) / matrixprofiler::std(query_window)
  sumwy <- sum(weight * query_window)
  sumwy2 <- sum(weight * (query_window^2))

  query_window <- rev(query_window)
  query_window[(window_size + 1):data_size] <- 0
  weight <- rev(weight)
  weight[(window_size + 1):data_size] <- 0

  queryw_fft <- stats::fft(weight * query_window)
  data_fft_queryw_fft <- data_fft * queryw_fft

  data_queryw <- stats::fft(data_fft_queryw_fft, inverse = TRUE) / length(data_fft_queryw_fft)

  last_product <- data_queryw[window_size:data_size]
  distance_profile <- data_pre - 2 * (last_product - sumwy * data_mean) / data_sd + sumwy2

  return(list(distance_profile = distance_profile, last_product = last_product))
}
