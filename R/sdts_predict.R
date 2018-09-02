#' Scalable Dictionary learning for Time Series (SDTS) prediction function
#'
#' @param model a model created by SDTS training function [sdts.train()].
#' @param data a `vector` of `numeric`. Time series.
#' @param window.size an `int`. The average sliding window size.
#'
#' @return Returns a `vector` of `logical` with predicted annotations.
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
#' tr_data <- test_data$train$data[subs]
#' tr_label <- test_data$train$label[subs]
#' te_data <- test_data$test$data[subs]
#' te_label <- test_data$test$label[subs]
#' model <- sdts.train(tr_data, tr_label, w, verbose = 0)
#' predict <- sdts.predict(model, te_data, round(mean(w)))
#' sdts.f.score(te_label, predict, 1)
#' \dontrun{
#' windows <- c(110, 220, 330)
#' model <- sdts.train(test_data$train$data, test_data$train$label, windows, verbose = 0)
#' predict <- sdts.predict(model, test_data$test$data, round(mean(windows)))
#' sdts.f.score(test_data$test$label, predict, 1)
#' }

sdts.predict <- function(model, data, window.size) {
  n.pat <- length(model$thold)
  anno.st <- list()
  data.size <- length(data)

  for (i in 1:n.pat) {
    pat.len <- length(model$pattern[[i]])
    pre <- mass.pre(data, data.size, window.size = pat.len)
    dist.pro <- mass(
      pre$data.fft, model$pattern[[i]], data.size, pat.len, pre$data.mean,
      pre$data.sd, mean(model$pattern[[i]]), std(model$pattern[[i]])
    )
    dist.pro <- Re(sqrt(dist.pro$distance.profile))
    anno <- dist.pro - model$thold[i]
    anno[anno >= 0] <- 0
    anno[anno < 0] <- -1
    anno <- -anno
    anno <- c(anno, rep(window.size, 0))
    anno.pad <- c(0, anno, 0)

    anno.st[[i]] <- which((anno.pad[1:(length(anno.pad) - 1)] - anno.pad[2:length(anno.pad)]) == -1) + 1
    anno.st[[i]] <- anno.st[[i]] - 1 ## ???
  }

  anno.st <- unlist(anno.st)
  anno.st <- sort(anno.st)

  i <- 1
  while (TRUE) {
    if (i >= length(anno.st)) {
      break
    }

    first.part <- anno.st[1:i]
    second.part <- anno.st[(i + 1):length(anno.st)]
    bad.st <- abs(second.part - anno.st[i]) < window.size
    second.part <- second.part[!bad.st]
    anno.st <- c(first.part, second.part)

    i <- i + 1
  }

  pred <- rep(FALSE, data.size - window.size + 1)
  anno.ed <- anno.st + window.size - 1

  for (i in 1:length(anno.st)) {
    pred[anno.st[i]:anno.ed[i]] <- TRUE
  }

  return(pred)
}

#' Computes the F-Score of a SDTS prediction
#'
#' Computes the F-Score of a SDTS prediction.
#'
#' `beta` is used to balance F-score towards recall (`>1`) or precision (`<1`).
#'
#' @param gtruth a `vector` of `logical`. Ground truth annotation.
#' @param pred a `vector` of `logical`. Predicted annotation from [sdts.predict()]
#' @param beta a `numeric`. See details. (default is `1`).
#'
#' @return Returns a `list` with `f.score`, `precision` and `recall`.
#'
#' @export
#'
#' @family SDTS
#' @references * Yeh C-CM, Kavantzas N, Keogh E. Matrix profile IV: Using Weakly Labeled Time Series
#'   to Predict Outcomes. Proc VLDB Endow. 2017 Aug 1;10(12):1802–12.
#' @references Website: <https://sites.google.com/view/weaklylabeled>
#' @examples
#' # This is a fast toy example and results are useless. For a complete result, run the code inside
#' #'Not run' section below.
#' w <- c(110, 220)
#' subs <- 11000:20000
#' tr_data <- test_data$train$data[subs]
#' tr_label <- test_data$train$label[subs]
#' te_data <- test_data$test$data[subs]
#' te_label <- test_data$test$label[subs]
#' model <- sdts.train(tr_data, tr_label, w, verbose = 0)
#' predict <- sdts.predict(model, te_data, round(mean(w)))
#' sdts.f.score(te_label, predict, 1)
#' \dontrun{
#' windows <- c(110, 220, 330)
#' model <- sdts.train(test_data$train$data, test_data$train$label, windows)
#' predict <- sdts.predict(model, test_data$test$data, round(mean(windows)))
#' sdts.f.score(test_data$test$label, predict, 1)
#' }
#'
sdts.f.score <- function(gtruth, pred, beta = 1) {
  if (length(pred) > length(gtruth)) {
    pred <- pred[1:length(gtruth)]
  } else if (length(pred) < length(gtruth)) {
    pred.tmp <- rep(FALSE, length(gtruth))
    pred.tmp[1:length(pred)] <- pred
    pred <- pred.tmp
  }

  pred.pad <- c(0, pred, 0)
  pred.st <- which((pred.pad[1:(length(pred.pad) - 1)] - pred.pad[2:length(pred.pad)]) == -1) + 1
  pred.ed <- which((pred.pad[1:(length(pred.pad) - 1)] - pred.pad[2:length(pred.pad)]) == 1)
  pred.st <- pred.st - 1
  pred.ed <- pred.ed - 1
  sub.len <- mode(pred.ed - pred.st + 1)

  is.tp <- rep(FALSE, length(pred.st))
  for (i in 1:length(pred.st)) {
    if (pred.ed[i] > length(gtruth)) {
      pred.ed[i] <- length(gtruth)
    }
    if (sum(gtruth[pred.st[i]:pred.ed[i]]) > 0.8 * sub.len) {
      is.tp[i] <- TRUE
    }
  }
  tp.pre <- sum(is.tp)

  gtruth.pad <- c(0, gtruth, 0)
  gtruth.st <- which((gtruth.pad[1:(length(gtruth.pad) - 1)] - gtruth.pad[2:length(gtruth.pad)]) == -1) + 1
  gtruth.ed <- which((gtruth.pad[1:(length(gtruth.pad) - 1)] - gtruth.pad[2:length(gtruth.pad)]) == 1)
  gtruth.st <- gtruth.st - 1
  gtruth.ed <- gtruth.ed - 1

  is.tp <- rep(FALSE, length(gtruth.st))
  for (i in 1:length(gtruth.st)) {
    if (gtruth.ed[i] > length(pred)) {
      gtruth.ed[i] <- length(gtruth)
    }
    if (sum(pred[gtruth.st[i]:gtruth.ed[i]]) > 0.8 * sub.len) {
      is.tp[i] <- TRUE
    }
  }
  tp.rec <- sum(is.tp)

  pre <- tp.pre / length(pred.st)
  rec <- tp.rec / length(gtruth.st)

  f.score <- (1 + beta^2) * (pre * rec) / ((beta^2) * pre + rec)

  return(list(f.score = f.score, precision = pre, recall = rec))
}

#' Calculates the mode of a vector
#'
#' @param x
#'
#' @return the mode
#' @keywords internal
#' @noRd

mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
