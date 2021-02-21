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
#' @param weight a `vector` of `numeric` or `NULL` with the same length of the `window_size`. This is
#' a MASS extension to weight the query.
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
dist_profile <- function(data, window_size, query = data, index = 1, params = NULL,
                         type = c("normalized", "non_normalized", "absolute", "weighted"), weights = NULL) {

  # Check time series classes -----------------------
  data <- convert_data(data)
  query <- convert_data(query)

  # Parse arguments ---------------------------------
  "!!!DEBUG Parsing data"
  checkmate::qassert(data, "N+")
  window_size <- checkmate::qassert(window_size, "x+[4,)")
  checkmate::qassert(query, "n+")
  checkmate::qassert(index, glue::glue("X1[1,{length(data) - window_size + 1}]"))
  "!!!DEBUG Parsing type"
  type <- match.arg(type)
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
      return(c(result, list(params = params)))
    },
    TRUE
  )

  if (is.null(query)) {
    "!!!DEBUG query is null"
    query <- data
  }

  # Computation ------------------------------------
  "!DEBUG Computation"
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
      result <- matrixprofiler::mass(params, data, query, index)
    },
    error = print
  )

  checkmate::qassert(params, "L+")
  checkmate::qassert(result, "L+")
}
