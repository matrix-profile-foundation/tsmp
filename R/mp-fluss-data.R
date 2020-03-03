#' Original data used in the FLUSS paper
#'
#' Contains two datasets used in FLUSS paper (1), first is TiltABP from (2), and second is
#' WalkJogRun from PAMAP's dataset (3)
#'
#' @docType data
#' @format A list containing:
#' \describe{
#'   \item{data}{one column matrix with the dataset's data}
#'   \item{gtruth}{a vector with the ground truth of semantic change according to provided dataset}
#'   \item{window}{window size used in original paper}
#' }
#' @source \url{https://sites.google.com/site/onlinesemanticsegmentation/}
#' @source \url{http://www.cs.ucr.edu/~eamonn/time_series_data/}
#'
#' @references * Gharghabi S, Ding Y, Yeh C-CM, Kamgar K, Ulanova L, Keogh E. Matrix Profile VIII:
#'   Domain Agnostic Online Semantic Segmentation at Superhuman Performance Levels. In: 2017 IEEE
#'   International Conference on Data Mining (ICDM). IEEE; 2017. p. 117-26.
#' @references * Heldt, T., Oefinger, M.B., Hoshiyama, M. and Mark, R.G., 2003, September.
#'   Circulatory response to passive and active changes in posture. In IEEE Computers in Cardiology,
#'   2003 (pp. 263-266).
#' @references * Reiss, A. and Stricker, D., 2012. Introducing a new benchmarked dataset for
#'   activity monitoring. In 16th International Symposium on Wearable Computers (ISWC), 2012, pages
#'   108-109. IEEE, 2012.
#' @keywords datasets
"mp_fluss_data"
