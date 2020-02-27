#' Write a TSMP object to JSON file.
#'
#' @param x a `MatrixProfile` or `PMP` object. If not, the `base::write()` function will be called.
#' @param file a `character` string with the output filename.
#' @param ... other arguments to be passed forward.
#'
#' @name write
#' @export
#' @examples
#'
#' result <- compute(mp_toy_data$data[, 1], 80)
#' \dontrun{
#' write(result, file = "output.json")
#' }
write <- function(x, ...) {
  UseMethod("write")
}

#' @export
write.default <- function(x, ...) {
  # little trick to allow using the base::write() when no class is matched
  pos <- which("package:tsmp" == search())
  get("write", pos = pos + 1, mode = "function", inherits = TRUE)(x, ...)
}

#' @name write
#' @export

write.MatrixProfile <- function(x, file, ...) {

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

  dgtz <- getOption("digits", 5)
  options(digits = 19)
  write(RJSONIO::toJSON(x,
    .inf = "Infinity", # default "" Infinity"
    .na = "NaN", # default "null"
    collapse = "", # default "\n"
    .withNames = TRUE, # default length(x) > 0 && length(names(x)) > 0
    asIs = NA # default NA
  ), file = file)
  options(digits = dgtz)
}

#' @name write
#' @export
write.PMP <- function(x, file, ...) {

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

  dgtz <- getOption("digits", 5)
  options(digits = 19)
  write(RJSONIO::toJSON(x,
    .inf = "Infinity", # default "" Infinity"
    .na = "NaN", # default "null"
    collapse = "", # default "\n"
    .withNames = TRUE, # default length(x) > 0 && length(names(x)) > 0
    asIs = NA # default NA
  ), file = file)
  options(digits = dgtz)
}

#' Read TSMP object from JSON file.
#'
#' @param x a `character` string with the input filename.
#' @param ... other arguments to be passed forward.
#'
#' @name read
#' @export
#'
#' @examples
#'
#' \dontrun{
#' result <- read("input.json")
#' }
read <- function(x, ...) {
  UseMethod("read")
}

#' @export
read.default <- function(x, ...) {
  mp <- RJSONIO::fromJSON(x, asText = FALSE, simplify = TRUE)

  if (mp$class == "MatrixProfile") {
    mp$mp <- as.matrix(mp$mp)
    mp$pi <- as.matrix(mp$pi + 1L)

    mp$data$ts <- as.vector(mp$data$ts)

    attributes(mp) <- list(
      names = names(mp),
      class = mp$class,
      join = mp$join,
      metric = mp$metric,
      algorithm = mp$algorithm
    )

    mp$metric <- NULL
    mp$join <- NULL
    mp$class <- NULL
    mp$algorithm <- NULL
  } else if (mp$class == "PMP") {
    mp$pmpi <- lapply(mp$pmpi, function(x) x + 1L)

    mp$data$ts <- as.vector(mp$data$ts)

    attributes(mp) <- list(
      names = names(mp),
      class = mp$class,
      join = mp$join,
      metric = mp$metric,
      algorithm = mp$algorithm
    )

    mp$metric <- NULL
    mp$join <- NULL
    mp$class <- NULL
    mp$algorithm <- NULL
  }

  return(mp)
}
