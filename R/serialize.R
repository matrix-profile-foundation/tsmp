ensure_rjsonio <- function() {
  if (!requireNamespace("RJSONIO", quietly = TRUE)) {
    stop("Package 'RJSONIO' is required for JSON serialization. Please install it.", call. = FALSE)
  }
}

#' Write a TSMP object to JSON file.
#'
#' @param x a `MatrixProfile` or `PMP` object.
#' @param file a `character` string with the output filename.
#' @param ... other arguments to be passed forward.
#'
#' @name write_json
#' @export
#' @examples
#'
#' result <- compute(mp_toy_data$data[, 1], 80)
#' \donttest{
#' write_json(result, file = file.path(tempdir(), "output.json"))
#' }

write_json <- function(x, ...) {
  UseMethod("write_json")
}

#' @export
write_json.default <- function(x, ...) {
  stop("JSON serialization requires a MatrixProfile or PMP object.", call. = FALSE)
}

#' @name write_json
#' @export
write_json.MatrixProfile <- function(x, file, ...) {
  ensure_rjsonio()

  # Parse arguments ---------------------------------
  checkmate::qassert(file, "S+")

  rplc <- function(y) {
    if (length(y) > 0) {
      return(I(y))
    } else {
      return(NULL)
    }
  }

  x$mp <- as.vector(x$mp)
  x$pi <- as.vector(x$pi) - 1L

  if (!is.null(x$motif)) {
    x$motif <- rapply(x$motif, rplc, how = "list")
  }

  if (!is.null(x$discord)) {
    x$discord <- rapply(x$discord, rplc, how = "list")
  }

  x$metric <- attr(x, "metric", TRUE)
  x$join <- attr(x, "join", TRUE)
  x$class <- attr(x, "class", TRUE)
  x$algorithm <- attr(x, "algorithm", TRUE)
  x$partial <- attr(x, "partial", exact = TRUE)
  x$origin <- attr(x, "origin", exact = TRUE)

  dgtz <- getOption("digits", 5)
  on.exit(options(digits = dgtz), add = TRUE)
  options(digits = 19)
  base::write(RJSONIO::toJSON(x,
    .inf = "Infinity", # default "" Infinity"
    .na = "NaN", # default "null"
    collapse = "", # default "\n"
    .withNames = TRUE, # default length(x) > 0 && length(names(x)) > 0
    asIs = NA # default NA
  ), file = file)
}

#' @name write_json
#' @export
write_json.PMP <- function(x, file, ...) {
  ensure_rjsonio()

  # Parse arguments ---------------------------------
  checkmate::qassert(file, "S+")

  rplc <- function(y) {
    if (length(y) > 0) {
      return(I(y))
    } else {
      return(NULL)
    }
  }

  x$pmpi <- lapply(x$pmpi, function(x) x - 1L)

  if (!is.null(x$motif)) {
    x$motif <- rapply(x$motif, rplc, how = "list")
  }

  if (!is.null(x$discord)) {
    x$discord <- rapply(x$discord, rplc, how = "list")
  }

  x$metric <- attr(x, "metric", TRUE)
  x$join <- attr(x, "join", TRUE)
  x$class <- attr(x, "class", TRUE)
  x$algorithm <- attr(x, "algorithm", TRUE)
  x$partial <- attr(x, "partial", exact = TRUE)
  x$origin <- attr(x, "origin", exact = TRUE)

  dgtz <- getOption("digits", 5)
  on.exit(options(digits = dgtz), add = TRUE)
  options(digits = 19)
  base::write(RJSONIO::toJSON(x,
    .inf = "Infinity", # default "" Infinity"
    .na = "NaN", # default "null"
    collapse = "", # default "\n"
    .withNames = TRUE, # default length(x) > 0 && length(names(x)) > 0
    asIs = NA # default NA
  ), file = file)
}

#' Read TSMP object from JSON file.
#'
#' @param x a `character` string with the input filename.
#' @param ... other arguments to be passed forward.
#'
#' @name read_json
#' @export
#'
#' @examples
#'
#' result <- compute(mp_toy_data$data[, 1], 80)
#' \donttest{
#' tempfile <- file.path(tempdir(), "output.json")
#' write_json(result, file = tempfile)
#' result <- read_json(tempfile)
#' }
read_json <- function(x, ...) {
  UseMethod("read_json")
}

#' @export
read_json.default <- function(x, ...) {
  ensure_rjsonio()
  mp <- RJSONIO::fromJSON(x, asText = FALSE, simplify = RJSONIO::Strict)

  if (mp$class == "MatrixProfile") {
    mp$mp <- as.vector(mp$mp)
    mp$pi <- as.vector(mp$pi + 1L)

    mp$data$ts <- as.vector(mp$data$ts)

    attributes(mp) <- list(
      names = names(mp),
      class = mp$class,
      join = mp$join,
      metric = mp$metric,
      algorithm = mp$algorithm,
      partial = mp$partial,
      origin = mp$origin
    )

    mp$metric <- NULL
    mp$join <- NULL
    mp$class <- NULL
    mp$algorithm <- NULL
    mp$partial <- NULL
    mp$origin <- NULL
  } else if (mp$class == "PMP") {
    mp$pmpi <- lapply(mp$pmpi, function(x) x + 1L)

    mp$data$ts <- as.vector(mp$data$ts)

    attributes(mp) <- list(
      names = names(mp),
      class = mp$class,
      join = mp$join,
      metric = mp$metric,
      algorithm = mp$algorithm,
      partial = mp$partial,
      origin = mp$origin
    )

    mp$metric <- NULL
    mp$join <- NULL
    mp$class <- NULL
    mp$algorithm <- NULL
    mp$partial <- NULL
    mp$origin <- NULL
  }

  return(mp)
}
