library(testthat)
if (!testthat:::on_cran()) {
  library(tsmp)
  Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")
  test_check("tsmp")
}
