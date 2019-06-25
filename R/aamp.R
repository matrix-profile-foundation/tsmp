#' matrix profile using the Euclidean distance
#'
#' @param data
#' @param window_size
#'
#' @return
#'
#' @examples
aamp <- function(data, window_size) {
  data_size <- length(data)
  s <- data_size - window_size
  matrix_profile <- rep(Inf, s)
  profile_index <- rep(-1, s)

  for (k in seq_len(s - 1)) {
    dist <- sum((data[1:window_size] - data[(k + 1):(k + window_size)])^2)

    if (dist < matrix_profile[1]) {
      matrix_profile[1] <- dist
      profile_index[1] <- k
    }

    if (dist < matrix_profile[k]) {
      matrix_profile[k] <- dist
      profile_index[k] <- 1
    }

    for (i in seq_len(s - k)) {
      kplusi <- k + i

      dist <- dist - (data[i] - data[kplusi])^2 + (data[window_size + i] - data[window_size + kplusi])^2

      if (matrix_profile[i] > dist) {
        matrix_profile[i] <- dist
        profile_index[i] <- kplusi
      }

      if (matrix_profile[kplusi] > dist) {
        matrix_profile[kplusi] <- dist
        profile_index[kplusi] <- i
      }
    }
  }

  matrix_profile <- sqrt(matrix_profile)
  return(list(mp = matrix_profile, pi = profile_index))
}

#' Z-Normalized Euclidean Distance
#'
#' @param data
#' @param window_size
#'
#' @return
#'
#' @examples
acamp <- function(data, window_size) {
  data_size <- length(data)
  s <- data_size - window_size
  matrix_profile <- rep(Inf, s)
  profile_index <- rep(-1, s)
  mm <- 1 / window_size

  sum_a <- sum(data[1:window_size]) # sum of the values in T1,window_size
  sum_a2 <- sum(data[1:window_size]^2) # sum of squares of the values in T1,window_size
  x1 <- data[1]
  xm <- data[1 + window_size]
  sum_b <- sum_a - x1 + xm
  sum_b2 <- sum_a2 - x1^2 + xm^2

  for (k in seq_len(s - 1)) {
    sum_ka <- sum_a
    sum_kb <- sum_b
    sum_ka2 <- sum_a2
    sum_kb2 <- sum_b2
    kplus1 <- 1 + k
    prod_c <- sum(data[kplus1:(k + window_size)] * data[1:window_size])
    z_product <- sum_ka * sum_kb - window_size * prod_c
    dist <- abs(z_product) * (z_product) / ((sum_ka2 - sum_ka^2 * mm) * (sum_kb2 - sum_kb^2 * mm))

    if (dist < matrix_profile[1]) {
      matrix_profile[1] <- dist
      profile_index[1] <- k
    }

    if (dist < matrix_profile[k]) {
      matrix_profile[k] <- dist
      profile_index[k] <- 1
    }


    sum_ka <- sum_ka - x1 + xm
    sum_ka2 <- sum_ka2 - x1^2 + xm^2
    xk <- data[kplus1]
    xkm <- data[window_size + kplus1]
    sum_kb <- sum_kb - xk + xkm
    sum_kb2 <- sum_kb2 - xk^2 + xkm^2
    sum_b <- sum_kb
    sum_b2 <- sum_kb2
    prod_c <- prod_c - x1 * xk + xm * xkm
    z_product <- sum_ka * sum_kb - window_size * prod_c
    dist <- abs(z_product) * (z_product) / ((sum_ka2 - sum_ka^2 * mm) * (sum_kb2 - sum_kb^2 * mm))

    if (matrix_profile[1] > dist) {
      profile_index[1] <- kplus1
      matrix_profile[1] <- dist
    }


    if (matrix_profile[kplus1] > dist) {
      profile_index[kplus1] <- 1
      matrix_profile[kplus1] <- dist
    }


    if ((s - k) > 1) {
      idxs <- 2:(s - k)
    } else {
      idxs <- 1
    }

    for (i in idxs) {
      xi <- data[i]
      xmi <- data[window_size + i]
      sum_ka <- sum_ka - xi + xmi
      sum_ka2 <- sum_ka2 - xi^2 + xmi^2
      iplusk <- i + k
      xik <- data[iplusk]
      xmik <- data[window_size + iplusk]
      sum_kb <- sum_kb - xik + xmik
      sum_kb2 <- sum_kb2 - xik^2 + xmik^2
      prod_c <- prod_c - xi * xik + xmi * xmik
      z_product <- sum_ka * sum_kb - window_size * prod_c
      dist <- abs(z_product) * (z_product) / ((sum_ka2 - sum_ka^2 * mm) * (sum_kb2 - sum_kb^2 * mm))

      if (matrix_profile[i] > dist) {
        profile_index[i] <- iplusk
        matrix_profile[i] <- dist
      }

      if (matrix_profile[iplusk] > dist) {
        profile_index[iplusk] <- i
        matrix_profile[iplusk] <- dist
      }
    }
  }

  mp_sign <- sign(matrix_profile)
  matrix_profile <- sqrt(2 * window_size + 2 * mp_sign * sqrt(mp_sign * matrix_profile))
  return(list(mp = matrix_profile, pi = profile_index))
}
