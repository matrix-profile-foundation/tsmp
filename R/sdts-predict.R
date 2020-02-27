#' Framework for Scalable Dictionary learning for Time Series (SDTS) prediction function
#'
#' This function trains a model that uses a dictionary to predict state changes. Differently from
#' [fluss()], it doesn't look for semantic changes (that may be several), but for binary states like
#' "on" or "off". Think for example that a human annotator is pressing a switch any time he thinks
#' that the recorded data is relevant, and releases the switch when he thinks the data is noise. This
#' algorithm will learn the switching points (even better) and try to predict using new data.
#'
#' @param model a model created by SDTS training function [sdts_train()].
#' @param data a `vector` of `numeric`. Time series.
#' @param window_size an `int`. The average sliding window size.
#'
#' @return Returns a `vector` of `logical` with predicted annotations.
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
#' model <- sdts_train(mp_test_data$train$data, mp_test_data$train$label, windows, verbose = 0)
#' predict <- sdts_predict(model, mp_test_data$test$data, round(mean(windows)))
#' sdts_score(predict, mp_test_data$test$label, 1)
#' }
#'
sdts_predict <- function(model, data, window_size) {
  n_pat <- length(model$thold)
  anno_st <- list()
  data_size <- length(data)

  for (i in 1:n_pat) {
    nn <- dist_profile(data, model$pattern[[i]])
    dist_pro <- sqrt(nn$distance_profile)
    anno <- dist_pro - model$thold[i]
    anno[anno >= 0] <- 0
    anno[anno < 0] <- -1
    anno <- -anno
    anno <- c(anno, rep(window_size, 0))
    anno_pad <- c(0, anno, 0)

    anno_st[[i]] <- which((anno_pad[1:(length(anno_pad) - 1)] - anno_pad[2:length(anno_pad)]) == -1)
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

  pred <- rep(FALSE, data_size - window_size + 1)
  anno_ed <- anno_st + window_size - 1

  for (i in seq_len(length(anno_st))) {
    pred[anno_st[i]:anno_ed[i]] <- TRUE
  }

  pred <- pred[seq_len(data_size - window_size + 1)]

  return(pred)
}

#' Computes the F-Score of a SDTS prediction
#'
#' Computes the F-Score of a SDTS prediction.
#'
#' `beta` is used to balance F-score towards recall (`>1`) or precision (`<1`).
#'
#' @param pred a `vector` of `logical`. Predicted annotation from [sdts_predict()]
#' @param gtruth a `vector` of `logical`. Ground truth annotation.
#' @param beta a `numeric`. See details. (default is `1`).
#'
#' @return Returns a `list` with `f_score`, `precision` and `recall`.
#'
#' @export
#'
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
#'
sdts_score <- function(pred, gtruth, beta = 1) {
  if (length(pred) > length(gtruth)) {
    pred <- pred[seq_len(length(gtruth))]
  } else if (length(pred) < length(gtruth)) {
    pred_tmp <- rep(FALSE, length(gtruth))
    pred_tmp[seq_len(length(pred))] <- pred
    pred <- pred_tmp
  }

  if (anyNA(gtruth)) {
    stop("`gtruth` contains NA values.", call. = FALSE)
  }

  if (anyNA(pred)) {
    stop("`pred` contains NA values.", call. = FALSE)
  }

  pred_pad <- c(0, pred, 0)
  pred_st <- which((pred_pad[1:(length(pred_pad) - 1)] - pred_pad[2:length(pred_pad)]) == -1) + 1
  pred_ed <- which((pred_pad[1:(length(pred_pad) - 1)] - pred_pad[2:length(pred_pad)]) == 1)

  pred_len <- min(length(pred_st), length(pred_ed))
  pred_st <- pred_st[seq_len(pred_len)]
  pred_ed <- pred_ed[seq_len(pred_len)]

  pred_st <- pred_st - 1
  pred_ed <- pred_ed - 1
  sub_len <- mode(pred_ed - pred_st + 1)

  is_tp <- rep(FALSE, length(pred_st))
  for (i in seq_len(length(pred_st))) {
    if (pred_ed[i] > length(gtruth)) {
      pred_ed[i] <- length(gtruth)
    }
    if (sum(gtruth[pred_st[i]:pred_ed[i]]) > 0.8 * sub_len) {
      is_tp[i] <- TRUE
    }
  }
  tp_pre <- sum(is_tp)

  gtruth_pad <- c(0, gtruth, 0)
  gtruth_st <- which((gtruth_pad[1:(length(gtruth_pad) - 1)] - gtruth_pad[2:length(gtruth_pad)]) == -1) + 1
  gtruth_ed <- which((gtruth_pad[1:(length(gtruth_pad) - 1)] - gtruth_pad[2:length(gtruth_pad)]) == 1)
  gtruth_st <- gtruth_st - 1
  gtruth_ed <- gtruth_ed - 1

  is_tp <- rep(FALSE, length(gtruth_st))
  for (i in seq_len(length(gtruth_st))) {
    if (gtruth_ed[i] > length(pred)) {
      gtruth_ed[i] <- length(gtruth)
    }
    if (sum(pred[gtruth_st[i]:gtruth_ed[i]]) > 0.8 * sub_len) {
      is_tp[i] <- TRUE
    }
  }
  tp_rec <- sum(is_tp)

  pre <- tp_pre / length(pred_st)
  rec <- tp_rec / length(gtruth_st)

  f_score <- (1 + beta^2) * (pre * rec) / ((beta^2) * pre + rec)

  return(invisible(list(f_score = f_score, precision = pre, recall = rec)))
}
