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
#' @param \dots Precomputed values from the first iteration. If not supplied, these values will be computed.
#' @param window_size an `int` or `NULL`. Sliding window size. See details.
#' @param method method that will be used to calculate the distance profile. See details.
#' @param index an `int`. Index of query window. See details.
#' @param k an `int` or `NULL`. Default is `NULL`. Defines the size of batch for MASS V3. Prefer to
#' use a power of 2. If `NULL`, it will be set automatically.
#' @param weight a `vector` of `numeric` or `NULL` with the same length of the `window_size`. This is
#' a MASS extension to weight the query.
#' @param paa a `numeric`. Default is `1`. Factor of PAA reduction (2 == half of size). This is a
#' MASS extension.
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
dist_profile <- function(data, query, ..., window_size = NULL, method = "v3", index = 1, k = NULL,
                         weight = NULL, paa = 1) {

  ## ---- Verify if method exists ----
  # set as v3 if no method is entered
  if (!is.na(pmatch(method, "v3"))) {
    method <- "v3"
  }
  valid_methods <- c(
    "v2", "v3", "weighted"
  )
  method <- pmatch(method, valid_methods)

  method <- paste0("mass_", valid_methods[method])

  if (is.na(method)) {
    stop("invalid similarity search method")
  }
  if (method == -1) {
    stop("ambiguous similarity search method")
  }

  # Grab parameters
  if (!missing(...) && is.list(...)) {
    params <- c(...)
  } else {
    params <- NULL
  }

  data <- as.vector(data)
  query <- as.vector(query)

  if (anyNA(query)) {
    ## ---- Query with Gap ----
    q1 <- NULL
    q2 <- NULL

    min_idx <- min(which(is.na(query))) - 1
    max_idx <- max(which(is.na(query))) + 1

    if (min_idx >= 4L) {
      q1 <- query[1:min_idx]
    }

    if (max_idx <= (window_size - 4)) {
      q2 <- query[max_idx:window_size]
    }

    if (anyNA(q1) || anyNA(q2)) {
      stop("Querying with gap only supports one gap.")
    }

    q1_len <- length(q1)
    q2_len <- length(q2)

    if (q1_len != q2_len) {
      warning("Warning: Result may be inconsistent if the size of the queries are different.")
    }

    pre1 <- mass_pre(data, window_size = q1_len)
    pre2 <- mass_pre(data, window_size = q2_len)

    result1 <- mass_v3(q1, data, pre1$window_size, pre1$data_size, pre1$data_mean, pre1$data_sd, mean(q1), std(q1))
    result2 <- mass_v3(q2, data, pre2$window_size, pre2$data_size, pre2$data_mean, pre2$data_sd, mean(q2), std(q2))
    result1 <- sqrt(result1$distance_profile)
    result2 <- sqrt(result2$distance_profile)

    pad <- rep(Inf, max_idx - 1)
    result2 <- c(pad, result2)
    result1 <- c(result1, rep(Inf, length(result2) - length(result1)))

    result <- list(distance_profile = (result1 + result2)^2)
  } else {
    ## ---- Non-GAP ----
    # set window_size
    window_size <- ifelse(is.null(window_size), length(query), window_size)

    if (paa > 1) {
      data <- paa(data, paa)
      query <- paa(query, paa)
      window_size <- floor(window_size / paa)
    }

    ## ---- First iteration with MASS ----
    if (is.null(params)) {
      params <- switch(method,
        mass_v2 = mass_pre(data, query, window_size),
        mass_v3 = c(mass_pre(data, query, window_size), list(data = data, k = k)),
        mass_weighted = mass_pre_w(data, query, window_size, weight)
      )
    } else {
      if (!is.null(params$par)) {
        params <- params$par
      }
      window_size <- params$window_size
    }

    pars <- params
    pars$query_mean <- pars$query_mean[index]
    pars$query_sd <- pars$query_sd[index]
    q_w <- query[index:(window_size + index - 1)]

    result <- do.call(method, c(list(q_w), pars))

    if (paa > 1) {
      result$distance_profile <- result$distance_profile * paa
      result$last_product <- result$last_product * paa
    }
  }

  return(c(result, list(par = params)))
}
