#' Scalable Dictionary learning for Time Series (SDTS) training function
#'
#' @details
#' `beta` is used to balance F-score towards recall (`>1`) or precision (`<1`). `verbose` changes
#' how much information is printed by this function; `0` means nothing, `1` means text, `2` means
#' text and sound.
#'
#' @param data a `vector` of `numeric`. Time series.
#' @param label a `vector` of `logical`. Annotations.
#' @param window_size an `int` or a `vector` of `int`. Sliding window sizes.
#' @param beta a `numeric`. See details. (default is `1`).
#' @param pat_max an `int`. Max number of shape features captured. (default is `Inf``).
#' @param parallel a `logical`. Use parallel computation inside (default is `TRUE`).
#' @param verbose an `int`. See details. (Default is `2`).
#'
#' @return Returns a list with the learned dictionary `score` (estimated score), `score_hist`
#'   (history of scores), `pattern` (shape features), `thold` (threshold values).
#'
#' @export
#' @family SDTS
#' @references * Yeh C-CM, Kavantzas N, Keogh E. Matrix profile IV: Using Weakly Labeled Time Series
#'   to Predict Outcomes. Proc VLDB Endow. 2017 Aug 1;10(12):1802–12.
#' @references Website: <https://sites.google.com/view/weaklylabeled>
#' @examples
#' # This is a fast toy example and results are useless. For a complete result, run the code inside
#' #'Not run' section below.
#' w <- c(110, 220)
#' subs <- 11000:20000
#' tr_data <- mp_test_data$train$data[subs]
#' tr_label <- mp_test_data$train$label[subs]
#' te_data <- mp_test_data$test$data[subs]
#' te_label <- mp_test_data$test$label[subs]
#' model <- sdts_train(tr_data, tr_label, w, verbose = 0)
#' predict <- sdts_predict(model, te_data, round(mean(w)))
#' sdts_score(te_label, predict, 1)
#' \dontrun{
#' windows <- c(110, 220, 330)
#' model <- sdts_train(mp_test_data$train$data, mp_test_data$train$label, windows)
#' predict <- sdts_predict(model, mp_test_data$test$data, round(mean(windows)))
#' sdts_score(mp_test_data$test$label, predict, 1)
#' }
sdts_train <- function(data, label, window_size, beta = 1, pat_max = Inf, parallel = TRUE, verbose = 2) {

  # transform data list into matrix
  if (is.matrix(data) || is.data.frame(data)) {
    if (is.data.frame(data)) {
      data <- as.matrix(data)
    } # just to be uniform
    if (ncol(data) > nrow(data)) {
      data <- t(data)
    }
    data_size <- nrow(data)
  } else if (is.vector(data)) {
    data_size <- length(data)
    # transform data into 1-col matrix
    data <- as.matrix(data) # just to be uniform
  } else {
    stop("Error: Unknown type of data. Must be: matrix, data.frame, vector or list.")
  }

  n_window_size <- length(window_size)

  # check input
  for (i in 1:n_window_size) {
    if (window_size[i] > (data_size / 2)) {
      stop("Error: Time series is too short relative to desired window size.")
    }
    if (window_size[i] < 4) {
      stop("Error: `window_size` must be at least 4.")
    }
  }

  # extract positive segment
  label_diff <- diff(c(0, label, 0))
  pos_st <- which(label_diff == 1) + 1
  pos_ed <- which(label_diff == -1)
  n_pos <- length(pos_st)
  pos_st <- pos_st - 1
  pos_ed <- pos_ed - 1
  pos <- list()

  for (i in 1:n_pos) {
    pos[[i * 2 - 1]] <- Inf
    pos[[i * 2]] <- data[pos_st[i]:pos_ed[i]]
  }

  pos <- unlist(pos)
  pos_alt_st <- which(is.infinite(pos)) + 1
  pos_alt_ed <- which(is.infinite(pos)) - 1
  pos_alt_ed <- c(pos_alt_ed[-1], length(pos))

  if (pos_alt_st[length(pos_alt_st)] > (length(pos) - min(window_size) + 1)) {
    pos_alt_st[length(pos_alt_st)] <- (length(pos) - min(window_size) + 1)
  }

  # run matrix profile on concatenated positive segment
  if (verbose > 0) {
    message("stage 1 of 3, compute matrix profile ...")
  }

  mat_pro <- list()

  for (i in 1:n_window_size) {
    if (parallel == TRUE) {
      mp <- stomp_par(pos, window_size = window_size[i], verbose = verbose)
    } else {
      mp <- stomp(pos, window_size = window_size[i], verbose = verbose)
    }
    mat_pro[[i]] <- mp$mp
  }

  # extract candidate
  candi <- list()
  candi_idx <- list()

  for (i in 1:n_window_size) {
    candi[[i]] <- list()
    candi_idx[[i]] <- rep(0, n_pos)
    candi_dist <- rep(0, n_pos)

    for (j in 1:n_pos) {
      temp <- mat_pro[[i]][pos_alt_st[j]:max(pos_alt_st[j], (pos_alt_ed[j] - window_size[i] + 1), na.rm = TRUE)]
      rlt_idx <- which.min(temp)
      candi_dist[j] <- temp[rlt_idx]

      alt_idx <- pos_alt_st[j] + rlt_idx - 1
      candi[[i]][[j]] <- pos[alt_idx:(alt_idx + window_size[i] - 1)]
      candi_idx[[i]][j] <- pos_st[j] + rlt_idx - 1
    }
    candi_dist <- sort(candi_dist, index.return = TRUE)
    candi[[i]] <- candi[[i]][candi_dist$ix]
    candi_idx[[i]] <- candi_idx[[i]][candi_dist$ix]
  }

  # evaluate each candidate
  candi_score <- list()
  candi_thold <- list()
  candi_pro <- list()
  candi_window_size <- list()
  tictac <- Sys.time()

  if (verbose > 0) {
    message("stage 2 of 3, evaluate individual candidate ...")
  }

  if (verbose > 0) {
    pb <- utils::txtProgressBar(min = 0, max = n_window_size * n_pos, style = 3, width = 80)
    on.exit(close(pb))
  }
  if (verbose > 1) {
    on.exit(beep(sounds[[1]]), TRUE)
  }

  for (i in 1:n_window_size) {
    candi_score[[i]] <- rep(0, n_pos)
    candi_thold[[i]] <- rep(0, n_pos)
    candi_pro[[i]] <- list()
    candi_window_size[[i]] <- rep(1, n_pos) * window_size[i]

    pre <- mass_pre(data, data_size, window_size = window_size[i])

    for (j in 1:n_pos) {
      dist_pro <- mass(pre$data_fft, candi[[i]][[j]], data_size, window_size[i], pre$data_mean, pre$data_sd, mean(candi[[i]][[j]]), std(candi[[i]][[j]]))
      dist_pro <- Re(sqrt(dist_pro$distance_profile))
      candi_pro[[i]][[j]] <- dist_pro
      exc_st <- max(1, candi_idx[[i]][j] - window_size[i])
      exc_ed <- min(length(dist_pro), candi_idx[[i]][j] + window_size[i])
      dist_pro[exc_st:exc_ed] <- Inf

      golden <- golden_section(dist_pro, label, pos_st, pos_ed, beta, window_size[i])

      candi_thold[[i]][j] <- golden$thold
      candi_score[[i]][j] <- golden$score

      if (verbose > 0) {
        utils::setTxtProgressBar(pb, ((i - 1) * n_pos + j))
      }
    }
  }

  tictac <- Sys.time() - tictac
  if (verbose > 0) {
    message(sprintf("\nFinished in %.2f %s", tictac, units(tictac)))
  }

  candi_pro_exp <- list()
  candi_exp <- list()

  for (i in 1:n_window_size) {
    candi_pro_exp[((i - 1) * n_pos + 1):(i * n_pos)] <- candi_pro[[i]]
    candi_exp[((i - 1) * n_pos + 1):(i * n_pos)] <- candi[[i]]
  }

  candi_pro <- candi_pro_exp
  candi <- candi_exp
  candi_score <- unlist(candi_score)
  candi_thold <- unlist(candi_thold)
  candi_idx <- unlist(candi_idx)
  candi_window_size <- unlist(candi_window_size)
  candi_score_sorted <- sort(signif(candi_score, 6), decreasing = TRUE, index.return = TRUE)
  candi_score <- candi_score_sorted$x
  order <- candi_score_sorted$ix

  candi_thold <- candi_thold[order]
  candi_idx <- candi_idx[order]
  candi_window_size <- candi_window_size[order]
  candi_pro <- candi_pro[order]
  candi <- candi[order]

  # check max pattern allowed
  pat_max <- min(pat_max, floor(n_pos * 0.5))
  if (pat_max < 2) {
    return(list(score = candi_score[1], score_hist = candi_score[1], pattern = list(candi[[1]]), thold = candi_thold[1]))
  }

  # check combined pattern
  max_window_size <- max(window_size)
  max_pro_len <- length(data) - min(window_size) + 1
  best_pat <- rep(FALSE, n_pos * n_window_size)
  best_score <- -Inf
  exc_mask <- rep(FALSE, max_pro_len)
  score_hist <- rep(Inf, n_pos * n_window_size)
  tictac <- Sys.time()

  if (verbose > 0) {
    message("stage 3 of 3, evaluate combination of candidates ...")
  }

  if (verbose > 0) {
    close(pb)
    pb <- utils::txtProgressBar(min = 0, max = pat_max * n_window_size * n_pos, style = 3, width = 80)
  }

  for (i in 1:pat_max) {
    pat_score <- rep(-Inf, n_pos * n_window_size)
    exc_mask_cur <- exc_mask
    exc_st <- rep(0, n_pos * n_window_size)
    exc_ed <- rep(0, n_pos * n_window_size)
    thold_cur <- list()

    for (j in 1:(n_pos * n_window_size)) {
      if (best_pat[j]) {
        next
      }

      best_pat_cur <- best_pat
      best_pat_cur[j] <- TRUE

      exc_st[j] <- max(1, candi_idx[j] - max_window_size)
      exc_ed[j] <- min(max_pro_len, candi_idx[j] + max_window_size)
      exc_mask_cur[exc_st[j]:exc_ed[j]] <- TRUE

      pro_cur <- candi_pro[best_pat_cur]
      pro_max <- -Inf
      pro_min <- Inf

      for (k in seq_len(length(pro_cur))) {
        pro_max <- max(max(pro_cur[[k]][!is.infinite(pro_cur[[k]])]), pro_max)
        pro_min <- min(min(pro_cur[[k]]), pro_min)
        pro_cur[[k]][exc_mask_cur] <- Inf
      }

      thold_cur[[j]] <- candi_thold[best_pat_cur]
      window_size_cur <- candi_window_size[best_pat_cur]

      iter <- 0
      score <- NULL
      while (TRUE) {
        iter <- iter + 1
        thold_old <- thold_cur[[j]]
        for (k in rev(seq_len(length(thold_cur[[j]])))) {
          gold <- golden_section_2(
            pro_cur,
            thold_cur[[j]],
            label,
            pos_st,
            pos_ed,
            beta,
            window_size_cur[k],
            k
          )
          thold_cur[[j]] <- gold$thold
          score <- gold$score
        }

        if ((iter > 200) || (mean(thold_cur[[j]] - thold_old) < ((pro_max - pro_min) * 0.001))) {
          break
        }
      }

      pat_score[j] <- score
      exc_mask_cur[exc_st[j]:exc_ed[j]] <- FALSE

      if (verbose > 0) {
        utils::setTxtProgressBar(pb, ((i - 1) * (n_pos * n_window_size) + j))
      }
    }
    if (verbose > 0) {
      utils::setTxtProgressBar(pb, ((i - 1) * (n_pos * n_window_size) + (n_pos * n_window_size)))
    }

    best_candi_idx <- which.max(pat_score)

    if ((pat_score[best_candi_idx] - best_score) > 0) {
      score_hist[i] <- pat_score[best_candi_idx]
      best_score <- pat_score[best_candi_idx]
      best_pat[best_candi_idx] <- TRUE
      candi_thold[best_pat] <- thold_cur[[best_candi_idx]]
      exc_mask[exc_st[best_candi_idx]:exc_ed[best_candi_idx]] <- TRUE
    } else {
      break
    }
  }

  if (verbose > 0) {
    utils::setTxtProgressBar(pb, (pat_max * n_pos * n_window_size))
  }

  tictac <- Sys.time() - tictac
  if (verbose > 0) {
    message(sprintf("\nFinished in %.2f %s", tictac, units(tictac)))
  }

  score_hist <- score_hist[!is.infinite(score_hist)]

  if (length(best_pat) == 1) {
    pattern <- list(candi[best_pat])
  } else {
    pattern <- candi[best_pat]
  }


  return(list(score = best_score, score_hist = score_hist, pattern = pattern, thold = candi_thold[best_pat]))
}

#' Computes the golden section for individual candidates
#'
#' @param dist_pro the candidate distance profile
#' @param label a vector with the data bool annotation
#' @param pos_st a vector with the starting points of label
#' @param pos_ed a vector with the ending points of label
#' @param beta a number that balance the F-Score. Beta > 1 towards recall, < towards precision
#' @param window_size an integer with the sliding window size
#'
#' @return Returns the best threshold and its F-Score
#'
#' @keywords internal
#' @noRd
#'
golden_section <- function(dist_pro, label, pos_st, pos_ed, beta, window_size) {
  golden_ratio <- (1 + sqrt(5)) / 2
  a_thold <- min(dist_pro)
  b_thold <- max(dist_pro[!is.infinite(dist_pro)])
  c_thold <- b_thold - (b_thold - a_thold) / golden_ratio
  d_thold <- a_thold + (b_thold - a_thold) / golden_ratio
  tol <- max((b_thold - a_thold) * 0.001, 0.0001)

  while (abs(c_thold - d_thold) > tol) {
    c_score <- compute_f_meas(label, pos_st, pos_ed, dist_pro, c_thold, window_size, beta)
    d_score <- compute_f_meas(label, pos_st, pos_ed, dist_pro, d_thold, window_size, beta)

    if (c_score$f_meas > d_score$f_meas) {
      b_thold <- d_thold
    } else {
      a_thold <- c_thold
    }

    c_thold <- b_thold - (b_thold - a_thold) / golden_ratio
    d_thold <- a_thold + (b_thold - a_thold) / golden_ratio
  }
  thold <- (a_thold + b_thold) * 0.5
  score <- compute_f_meas(label, pos_st, pos_ed, dist_pro, thold, window_size, beta)

  return(list(thold = thold, score = score$f_meas))
}

#' Computes the golden section for combined candidates

#' @param dist_pro the candidate distance profile
#'
#' @param thold a number with the threshold used to calculate the F-Score
#' @param label a vector with the data bool annotation
#' @param pos_st a vector with the starting points of label
#' @param pos_ed a vector with the ending points of label
#' @param beta a number that balance the F-Score. Beta > 1 towards recall, < towards precision
#' @param window_size an integer with the sliding window size
#' @param fit_idx an integer with the index of the current threshold
#'
#' @return Returns the best threshold and its F-Score
#'
#' @keywords internal
#' @noRd

golden_section_2 <- function(dist_pro, thold, label, pos_st, pos_ed, beta, window_size, fit_idx) {
  golden_ratio <- (1 + sqrt(5)) / 2
  a_thold <- min(dist_pro[[fit_idx]], na.rm = TRUE)
  b_thold <- max(dist_pro[[fit_idx]][!is.infinite(dist_pro[[fit_idx]])], na.rm = TRUE)
  c_thold <- b_thold - (b_thold - a_thold) / golden_ratio
  d_thold <- a_thold + (b_thold - a_thold) / golden_ratio
  tol <- max((b_thold - a_thold) * 0.001, 0.0001)

  while (abs(c_thold - d_thold) > tol) {
    c_thold_combined <- thold
    d_thold_combined <- thold
    c_thold_combined[fit_idx] <- c_thold
    d_thold_combined[fit_idx] <- d_thold

    c_score <- compute_f_meas(label, pos_st, pos_ed, dist_pro, c_thold_combined, window_size, beta)
    d_score <- compute_f_meas(label, pos_st, pos_ed, dist_pro, d_thold_combined, window_size, beta)

    if (c_score$f_meas > d_score$f_meas) {
      b_thold <- d_thold
    } else {
      a_thold <- c_thold
    }

    c_thold <- b_thold - (b_thold - a_thold) / golden_ratio
    d_thold <- a_thold + (b_thold - a_thold) / golden_ratio
  }
  thold[fit_idx] <- (a_thold + b_thold) * 0.5
  score <- compute_f_meas(label, pos_st, pos_ed, dist_pro, thold, window_size, beta)

  # beta = 2;   emphacise recall
  # beta = 0.5; emphacise precision
  return(list(thold = thold, score = score$f_meas))
}

#' Computes de F-Score
#'
#' @param label a vector with the data bool annotation
#' @param pos_st a vector with the starting points of label
#' @param pos_ed a vector with the ending points of label
#' @param dist_pro the distance profile of the data
#' @param thold a number with the threshold used to compute the prediction
#' @param window_size an integer with the sliding window size
#' @param beta a number that balance the F-Score. Beta > 1 towards recall, < towards precision
#'
#' @return Returns the F-Score, precision and recall values
#'
#' @keywords internal
#' @noRd

compute_f_meas <- function(label, pos_st, pos_ed, dist_pro, thold, window_size, beta) {
  # generate annotation curve for each pattern
  if (is.list(dist_pro)) {
    anno_st <- list()
    n_pat <- length(dist_pro)

    for (i in 1:n_pat) {
      annor <- dist_pro[[i]] - thold[i]
      annor[annor > 0] <- 0
      annor[annor < 0] <- -1
      annor <- -annor
      anno_st[[i]] <- which(diff(c(0, annor, 0)) == 1) + 1
      anno_st[[i]] <- anno_st[[i]] - 1
    }

    anno_st <- unlist(anno_st)
    anno_st <- sort(anno_st)

    i <- 1
    while (TRUE) {
      if (i >= length(anno_st)) {
        break
      }

      first_part <- anno_st[1:i]
      second_part <- anno_st[(i + 1):length(anno_st)]
      bad_st <- abs(second_part - anno_st[i]) < window_size

      second_part <- second_part[!bad_st]
      anno_st <- c(first_part, second_part)
      i <- i + 1
    }

    anno_ed <- anno_st + window_size - 1
  } else {
    anno <- dist_pro - thold
    anno[anno > 0] <- 0
    anno[anno < 0] <- -1
    anno <- -anno

    anno_st <- which(diff(c(0, anno, 0)) == 1) + 1
    anno_ed <- anno_st + window_size - 1
    anno_st <- anno_st - 1
    anno_ed <- anno_ed - 1
  }

  anno <- rep(FALSE, length(label))

  for (i in seq_len(length(anno_st))) {
    anno[anno_st[i]:anno_ed[i]] <- 1
  }

  is.tp <- rep(FALSE, length(anno_st))

  for (i in seq_len(length(anno_st))) {
    if (anno_ed[i] > length(label)) {
      anno_ed[i] <- length(label)
    }
    if (sum(label[anno_st[i]:anno_ed[i]]) > (0.8 * window_size)) {
      is.tp[i] <- TRUE
    }
  }
  tp_pre <- sum(is.tp)

  is.tp <- rep(FALSE, length(pos_st))
  for (i in seq_len(length(pos_st))) {
    if (sum(anno[pos_st[i]:pos_ed[i]]) > (0.8 * window_size)) {
      is.tp[i] <- TRUE
    }
  }
  tp_rec <- sum(is.tp)

  pre <- tp_pre / length(anno_st)
  rec <- tp_rec / length(pos_st)

  f_meas <- (1 + beta^2) * (pre * rec) / ((beta^2) * pre + rec)
  if (is.na(f_meas)) {
    f_meas <- 0
  }
  return(list(f_meas = f_meas, pre = pre, rec = rec))
}
