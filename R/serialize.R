#' Title
#'
#' @param mp
#' @param file
#'
#' @return
#' @export
#'
#' @examples
write_mp <- function(mp, file) {

  mp$mp <- as.vector(mp$mp)
  mp$pi <- as.vector(mp$pi)

  write(RJSONIO::toJSON(mp,
    .inf = "Infinity", # default "" Infinity"
    .na = "NaN", # default "null"
    collapse = "", # default "\n"
    .withNames = TRUE, # default length(x) > 0 && length(names(x)) > 0
    asIs = NA # default NA
  ), file = file)
}

read_mp <- function(file) {
  js <- RJSONIO::fromJSON(file, asText = FALSE, simplify = TRUE)

  return(js)
}
