library(testthat)
library(tsmp)

if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
  Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")
}

test_check("tsmp")
