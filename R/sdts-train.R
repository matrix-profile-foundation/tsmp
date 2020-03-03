#' Framework for Scalable Dictionary learning for Time Series (SDTS) training function
#'
#' This function trains a model that uses a dictionary to predict state changes. Differently from
#' [fluss()], it doesn't look for semantic changes (that may be several), but for binary states like
#' "on" or "off". Think for example that a human annotator is pressing a switch any time he thinks
#' that the recorded data is relevant, and releases the switch when he thinks the data is noise. This
#' algorithm will learn the switching points (even better) and try to predict using new data.
#'
#' @details
#' `beta` is used to balance F-score towards recall (`>1`) or precision (`<1`). `verbose` changes
#' how much information is printed by this function; `0` means nothing, `1` means text, `2` adds the
#' progress bar, `3` adds the finish sound.
#'
#' @param data a `vector` of `numeric`. Time series.
#' @param label a `vector` of `logical`. Annotations.
#' @param window_size an `int` or a `vector` of `int`. Sliding window sizes.
#' @param beta a `numeric`. See details. (default is `1`).
#' @param pat_max an `int`. Max number of shape features captured. (default is `Inf`).
#' @param parallel a `logical`. Use parallel computation inside (default is `FALSE`).
#' @param verbose an `int`. See details. (Default is `2`).
#'
#' @return Returns a list with the learned dictionary `score` (estimated score), `score_hist`
#'   (history of scores), `pattern` (shape features), `thold` (threshold values).
#'
#' @export
#' @family Scalable Dictionaries
#' @references * Yeh C-CM, Kavantzas N, Keogh E. Matrix profile IV: Using Weakly Labeled Time Series
#'   to Predict Outcomes. Proc VLDB Endow. 2017 Aug 1;10(12):1802-12.
#' @references Website: <https://sites.google.com/view/weaklylabeled>
#' @examples
#' # This is a fast toy example and results are useless. For a complete result, run the code inside
#' #' Not run' section below.
#' w <- c(110, 220)
#' subs <- 11000:20000
#' tr_data <- mp_test_data$train$data[subs]
#' tr_label <- mp_test_data$train$label[subs]
#' te_data <- mp_test_data$test$data[subs]
#' te_label <- mp_test_data$test$label[subs]
#' model <- sdts_train(tr_data, tr_label, w, verbose = 0)
#' predict <- sdts_predict(model, te_data, round(mean(w)))
#' sdts_score(predict, te_label, 1)
#' \dontrun{
#' windows <- c(110, 220, 330)
#' model <- sdts_train(mp_test_data$train$data, mp_test_data$train$label, windows)
#' predict <- sdts_predict(model, mp_test_data$test$data, round(mean(windows)))
#' sdts_score(predict, mp_test_data$test$label, 1)
#' }
sdts_train <- function(data, label, window_size, beta = 1, pat_max = Inf, parallel = FALSE, verbose = getOption("tsmp.verbose", 2)) {

  # transform data list into matrix ----
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
    stop("Unknown type of data. Must be: matrix, data.frame, vector or list.")
  }

  n_window_size <- length(window_size)

  # check input ----
  for (i in 1:n_window_size) {
    if (window_size[i] > (data_size / 2)) {
      stop("Time series is too short relative to desired window size.")
    }
    if (window_size[i] < 4) {
      stop("`window_size` must be at least 4.")
    }
  }

  # extract positive segment ----
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
  pos_alt_st <- list()
  for (i in 1:n_window_size) {
    pos_alt_st[[i]] <- which(is.infinite(pos)) + 1
    max_pos_idx <- pos_alt_st[[i]] > (length(pos) - window_size[i] + 1)
    pos_alt_st[[i]][max_pos_idx] <- (length(pos) - window_size[i] + 1)
  }

  pos_alt_ed <- which(is.infinite(pos)) - 1
  pos_alt_ed <- c(pos_alt_ed[-1], length(pos))



  # run matrix profile on concatenated positive segment ----
  if (verbose > 0) {
    message("Stage 1 of 3, compute matrix profile...")
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

  # extract candidate ----
  candi <- list()
  candi_idx <- list()

  for (i in 1:n_window_size) {
    candi[[i]] <- list()
    candi_idx[[i]] <- rep(0, n_pos)
    candi_dist <- rep(0, n_pos)

    for (j in 1:n_pos) {
      temp <- mat_pro[[i]][pos_alt_st[[i]][j]:max(pos_alt_st[[i]][j], (pos_alt_ed[j] - window_size[i] + 1), na.rm = TRUE)]
      rlt_idx <- which.min(temp)
      if (length(temp[rlt_idx]) == 0) {
        print("Zero")
      }

      candi_dist[j] <- temp[rlt_idx]

      alt_idx <- pos_alt_st[[i]][j] + rlt_idx - 1
      candi[[i]][[j]] <- pos[alt_idx:(alt_idx + window_size[i] - 1)]
      candi_idx[[i]][j] <- pos_st[j] + rlt_idx - 1
    }
    candi_dist <- sort(candi_dist, index.return = TRUE)
    candi[[i]] <- candi[[i]][candi_dist$ix]
    candi_idx[[i]] <- candi_idx[[i]][candi_dist$ix]
  }

  ### CHECK oK
  # evaluate each candidate ----
  candi_score <- list()
  candi_thold <- list()
  candi_pro <- list()
  candi_window_size <- list()
  tictac <- Sys.time()

  if (verbose > 0) {
    message("Stage 2 of 3, evaluate individual candidates...")
  }

  if (verbose > 1) {
    pb <- progress::progress_bar$new(
      format = "SDTS-Train [:bar] :percent at :tick_rate it/s, elapsed: :elapsed, eta: :eta",
      clear = FALSE, total = n_window_size * n_pos, width = 80
    )
  }
  if (verbose > 2) {
    on.exit(beep(sounds[[1]]), TRUE)
  }

  for (i in 1:n_window_size) {
    candi_score[[i]] <- rep(0, n_pos)
    candi_thold[[i]] <- rep(0, n_pos)
    candi_pro[[i]] <- list()
    candi_window_size[[i]] <- rep(1, n_pos) * window_size[i]

    for (j in 1:n_pos) {
      nn <- dist_profile(data, candi[[i]][[j]], window_size = window_size[i])
      dist_pro <- sqrt(nn$distance_profile)
      candi_pro[[i]][[j]] <- dist_pro
      exc_st <- max(1, candi_idx[[i]][j] - window_size[i])
      exc_ed <- min(length(dist_pro), candi_idx[[i]][j] + window_size[i])
      dist_pro[exc_st:exc_ed] <- Inf

      golden <- golden_section(dist_pro, label, pos_st, pos_ed, beta, window_size[i])

      candi_thold[[i]][j] <- golden$thold
      candi_score[[i]][j] <- golden$score

      if (verbose > 1) {
        pb$tick()
      }
    }
  }

  tictac <- Sys.time() - tictac
  if (verbose > 0) {
    message(sprintf("Finished in %.2f %s", tictac, units(tictac)))
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

  # check max pattern allowed ----
  pat_max <- min(pat_max, floor(n_pos * 0.5))
  if (pat_max < 2) {
    return(list(score = candi_score[1], score_hist = candi_score[1], pattern = list(candi[[1]]), thold = candi_thold[1]))
  }
  # CHECK OK
  # check combined pattern ----
  max_window_size <- max(window_size)
  max_pro_len <- length(data) - min(window_size) + 1
  best_pat <- rep(FALSE, n_pos * n_window_size)
  best_score <- -Inf
  exc_mask <- rep(FALSE, max_pro_len)
  score_hist <- rep(Inf, n_pos * n_window_size)
  tictac <- Sys.time()

  if (verbose > 0) {
    message("Stage 3 of 3, evaluate combination of candidates...")
  }

  if (verbose > 1) {
    pb <- progress::progress_bar$new(
      format = "SDTS-Train [:bar] :percent at :tick_rate it/s, elapsed: :elapsed, eta: :eta",
      clear = FALSE, total = pat_max * n_window_size * n_pos, width = 80
    )
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
        # exc_mask_cur is trimmed because if not, will insert NA in pro_cur
        pro_cur[[k]][exc_mask_cur[seq_len(length(pro_cur[[k]]))]] <- Inf
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

        if (iter > 200) {
          break
        }

        if (!anyNA(c(thold_cur[[j]], thold_old, pro_max, pro_min))) {
          if ((mean(thold_cur[[j]] - thold_old) < ((pro_max - pro_min) * 0.001))) {
            break
          }
        }
      }

      pat_score[j] <- score
      exc_mask_cur[exc_st[j]:exc_ed[j]] <- FALSE

      if (verbose > 1) {
        pb$tick()
      }
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

  if (verbose > 1) {
    pb$update(ratio = 1)
  }

  tictac <- Sys.time() - tictac
  if (verbose > 0) {
    message(sprintf("Finished in %.2f %s", tictac, units(tictac)))
  }

  score_hist <- score_hist[!is.infinite(score_hist)]

  if (length(best_pat) == 1) {
    pattern <- list(candi[best_pat])
  } else {
    pattern <- candi[best_pat]
  }


  return(list(score = best_score, score_hist = score_hist, pattern = pattern, thold = candi_thold[best_pat]))
}
