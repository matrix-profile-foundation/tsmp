#' Original data used in the Salient Subsequences demo
#'
#' This is the Meat dataset from UCR Archive modified for Salient discovery.
#' The original data is mixed with Random Walks and the algorithm must pick only the originals.
#'
#' @docType data
#' @format `original` is the original dataset with 60+60 observations mixed with 120 random walks:
#' \describe{
#'   \item{data}{240 time series with length of 448 each.}
#'   \item{labels}{label of each time series, `-666` means a random walk.}
#'   \item{sub_len}{size of sliding window.}
#' }
#'
#' `sub` is the original dataset embedded in random walks:
#'
#' \describe{
#'   \item{data}{One time series with length of 107520.}
#'   \item{labels}{label of each original data.}
#'   \item{labels_idx}{starting point where the original data was placed.}
#'   \item{sub_len}{size of sliding window.}
#' }
#'
#' @source \url{http://www.cs.ucr.edu/~eamonn/time_series_data/}
#'
#' @references * Yeh CCM, Van Herle H, Keogh E. Matrix profile III: The matrix profile allows
#'   visualization of salient subsequences in massive time series. Proc - IEEE Int Conf Data Mining,
#'   ICDM. 2017;579-88.
#' @references * Hu B, Rakthanmanon T, Hao Y, Evans S, Lonardi S, Keogh E. Discovering the Intrinsic
#'   Cardinality and Dimensionality of Time Series Using MDL. In: 2011 IEEE 11th International
#'   Conference on Data Mining. IEEE; 2011. p. 1086-91.
#' @references Website: <https://sites.google.com/site/salientsubs/>
#' @keywords datasets
"mp_meat_data"
