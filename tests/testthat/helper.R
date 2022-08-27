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

helper_basics <- function() {
  ref_data <- mp_toy_data$data[[1]]
  query_data <- mp_toy_data$data[[1]]
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
