#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
write <- function(x, ...) {
  UseMethod("write")
}

write.default <- function(...) {
  pos <- which("package:tsmp" == search())
  get("write", pos = pos + 1, mode = "function", inherits = TRUE)(...)
}

write.MatrixProfile <- function(x, ...) {

  pars <- list(...)

  rplc <- function(y) {
    if(length(y) > 0) {
      return(I(y))
    } else {
      return(NULL)
    }
  }

  x$mp <- as.vector(x$mp)
  x$pi <- as.vector(x$pi)
  x$motif <- rapply(x$motif, rplc, how = "list")
  x$discord <- rapply(x$discord, rplc, how = "list")

  dgtz <- getOption("digits", 5)
  options(digits = 19)
  write(RJSONIO::toJSON(x,
    .inf = "Infinity", # default "" Infinity"
    .na = "NaN", # default "null"
    collapse = "", # default "\n"
    .withNames = TRUE, # default length(x) > 0 && length(names(x)) > 0
    asIs = NA # default NA
  ), file = pars$file)
  options(digits = dgtz)
}

write.PanMatrixProfile <- function(x, ...) {


}

#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
read <- function(x, ...) {
  UseMethod("read")

}

read.default <- function(x, ...) {
  mp <- RJSONIO::fromJSON(x, asText = FALSE, simplify = TRUE)

  mp$mp <- as.matrix(mp$mp)
  mp$pi <- as.matrix(mp$pi)
  mp$data$ts <- as.vector(mp$data$ts)

  return(mp)
}
