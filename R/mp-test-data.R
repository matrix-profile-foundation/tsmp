#' Original data used in the STDS demo
#'
#' A synthetic dataset base on TRACE dataset and used as Stress Test to STDS algorithm. The TRACE
#' dataset used here is originally from (1), and the version distributed here is from (2)
#'
#' @docType data
#' @format A list of matrices with 215010 rows and 1 dimension:
#' \describe{
#'   \item{train$data}{training data}
#'   \item{train$label}{label for training data}
#'   \item{test$data}{test data}
#'   \item{test$label}{label for test data}
#' }
#' @source \url{https://sites.google.com/view/weaklylabeled}
#' @source \url{http://www.cs.ucr.edu/~eamonn/time_series_data/}
#'
#' @references * Roverso, D., Multivariate temporal classification by windowed wavelet decomposition
#'   and recurrent neural networks, in 3rd ANS Int'l Topical Meeting on Nuclear Plant
#'   Instrumentation, Control and Human-Machine Interface, vol. 20, Washington, DC, USA, 2000.
#' @references * Yeh C-CM, Kavantzas N, Keogh E. Matrix profile IV: Using Weakly Labeled Time Series
#'   to Predict Outcomes. Proc VLDB Endow. 2017 Aug 1;10(12):1802-12.
#' @keywords datasets
"mp_test_data"
