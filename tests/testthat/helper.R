new_rng_snapshots <- utils::compareVersion("3.6.0", as.character(getRversion())) > 0

helper_objects <- function() {
  ts <- mp_toy_data$data[[1]]
  query <- mp_toy_data$data[[2]]

  list(
    ts = ts,
    query = query
  )
}

helper_annotation <- function() {
  data <- mp_test_data$train$data[1:1000]

  list(
    data = data
  )
}

helper_dist_profile <- function() {
  ref_data <- mp_toy_data$data[[1]]
  query_data <- mp_toy_data$data[[2]]
  query_gap <- c(10:1, rep(NA, 10), 10:20)
  d_size <- length(ref_data)
  q_size <- length(query_data)
  w <- 30

  list(
    ref_data = ref_data,
    query_data = query_data,
    query_gap = query_gap,
    d_size = d_size,
    q_size = q_size,
    w = w
  )
}

helper_contrast <- function() {
  data1 <- mp_toy_data$data[[1]]
  data2 <- mp_toy_data$data[[2]]
  w <- 50

  list(
    data1 = data1,
    data2 = data2,
    w = w
  )
}


skip_if_not_ci <- function() {
  ci_providers <- c("GITHUB_ACTIONS", "TRAVIS", "APPVEYOR")
  ci <- any(toupper(Sys.getenv(ci_providers)) == "TRUE")
  if (ci) {
    return(invisible(TRUE))
  }
  skip("Not on GitHub Actions, Travis, or Appveyor")
}

skip_if_no_git_user <- function() {
  user_name <- git_cfg_get("user.name")
  user_email <- git_cfg_get("user.email")
  user_name_exists <- !is.null(user_name)
  user_email_exists <- !is.null(user_email)
  if (user_name_exists && user_email_exists) {
    return(invisible(TRUE))
  }
  skip("No Git user configured")
}

# CRAN's mac builder sets $HOME to a read-only ram disk, so tests can fail if
# you even tickle something that might try to lock its own config file during
# the operation (e.g. git) or if you simply test for writeability
skip_on_cran_macos <- function() {
  sysname <- tolower(Sys.info()[["sysname"]])
  on_cran <- !identical(Sys.getenv("NOT_CRAN"), "true")
  if (on_cran && sysname == "darwin") {
    skip("On CRAN and on macOS")
  }
  invisible(TRUE)
}

# with_mock <- function(..., .parent = parent.frame()) {
#   mockr::with_mock(..., .parent = .parent, .env = "usethis")
# }

expect_error_free <- function(...) {
  expect_error(..., regexp = NA)
}

is_build_ignored <- function(pattern, ..., base_path = usethis::proj_get()) {
  lines <- xfun::read_utf8(fs::path(base_path, ".Rbuildignore"))
  length(grep(pattern, x = lines, fixed = TRUE, ...)) > 0
}

# test_file <- function(fname) {
#   testthat::test_path("ref", fname)
# }

expect_proj_file <- function(...) {
  expect_true(fs::file_exists(usethis::proj_path(...)))
}
expect_proj_dir <- function(...) {
  expect_true(fs::dir_exists(usethis::proj_path(...)))
}
