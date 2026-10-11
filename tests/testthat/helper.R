expect_hash <- function(object, hash) {
  testthat::expect_equal(substr(digest::digest(object), 1, 10), hash)
}
