#' Scalable Dictionary learning for Time Series (SDTS) training function.
#'
#' Scalable Dictionary learning for Time Series (SDTS) training function.
#'
#' `beta` is used to balance F-score towards recall (`>1`) or precision (`<1`).
#'
#' @param data a `vector` of `numeric`. Time series.
#' @param label a `vector` of `logical`. Annotations.
#' @param window.size an `int` or a `vector` of `int`. Sliding window sizes.
#' @param beta a `numeric`. See details. (default is `1`).
#' @param pat.max an `int`. Max number of shape features captured. (default is `Inf``).
#' @param parallel a `logical`. Use parallel computation inside (default is `TRUE`).
#'
#' @return Returns a list with the learned dictionary
#'    `score` (estimated score), `score.hist` (history of scores),
#'    `pattern` (shape features), `thold` (threshold values).
#'
#' @export
#' @family SDTS
#'
#' @examples
#' \dontrun{
#' windows <- c(110, 220, 330)
#' model <- sdts.train(test_data$train$data, test_data$train$label, windows)
#' predict <- sdts.predict(model, test_data$test$data, round(mean(windows)))
#' sdts.f.score(test_data$test$label, predict, 1)
#' }
sdts.train <- function(data, label, window.size, beta = 1, pat.max = Inf, parallel = TRUE) {

  ## transform data list into matrix
  if (is.matrix(data) || is.data.frame(data)) {
    if (is.data.frame(data)) {
      data <- as.matrix(data)
    } # just to be uniform
    if (ncol(data) > nrow(data)) {
      data <- t(data)
    }
    data.size <- nrow(data)
  } else if (is.vector(data)) {
    data.size <- length(data)
    # transform data into 1-col matrix
    data <- as.matrix(data) # just to be uniform
  } else {
    stop("Unknown type of data. Must be: matrix, data.frame, vector or list")
  }

  n.window.size <- length(window.size)

  ## check input
  for (i in 1:n.window.size) {
    if (window.size[i] > (data.size / 2)) {
      stop("Error: Time series is too short relative to desired subsequence length")
    }
    if (window.size[i] < 4) {
      stop("Error: Subsequence length must be at least 4")
    }
  }

  ## extract positive segment
  label.diff <- diff(c(0, label, 0))
  pos.st <- which(label.diff == 1) + 1
  pos.ed <- which(label.diff == -1)
  n.pos <- length(pos.st)
  pos.st <- pos.st - 1
  pos.ed <- pos.ed - 1
  pos <- list()

  for (i in 1:n.pos) {
    pos[[i * 2 - 1]] <- Inf
    pos[[i * 2]] <- data[pos.st[i]:pos.ed[i]]
  }

  pos <- unlist(pos)
  pos.alt.st <- which(is.infinite(pos)) + 1
  pos.alt.ed <- which(is.infinite(pos)) - 1
  pos.alt.ed <- c(pos.alt.ed[-1], length(pos))

  if (pos.alt.st[length(pos.alt.st)] > (length(pos) - min(window.size) + 1)) {
    pos.alt.st[length(pos.alt.st)] <- (length(pos) - min(window.size) + 1)
  }

  ## run matrix profile on concatenated positive segment
  message("stage 1 of 3, compute matrix profile ...")

  mat.pro <- list()

  for (i in 1:n.window.size) {
    if (parallel == TRUE) {
      mp <- mstomp.par(pos, window.size[i])
    } else {
      mp <- mstomp(pos, window.size[i])
    }
    mat.pro[[i]] <- mp$mp
  }

  ## extract candidate
  candi <- list()
  candi.idx <- list()

  for (i in 1:n.window.size) {
    candi[[i]] <- list()
    candi.idx[[i]] <- rep(0, n.pos)
    candi.dist <- rep(0, n.pos)

    for (j in 1:n.pos) {
      temp <- mat.pro[[i]][pos.alt.st[j]:max(pos.alt.st[j], (pos.alt.ed[j] - window.size[i] + 1), na.rm = TRUE)]
      rlt.idx <- which.min(temp)
      candi.dist[j] <- temp[rlt.idx]

      alt.idx <- pos.alt.st[j] + rlt.idx - 1
      candi[[i]][[j]] <- pos[alt.idx:(alt.idx + window.size[i] - 1)]
      candi.idx[[i]][j] <- pos.st[j] + rlt.idx - 1
    }
    candi.dist <- sort(candi.dist, index.return = TRUE)
      # sort(signif(candi.dist, 6), index.return = TRUE)
    candi[[i]] <- candi[[i]][candi.dist$ix]
    candi.idx[[i]] <- candi.idx[[i]][candi.dist$ix]
  }

  ## evaluate each candidate
  candi.score <- list()
  candi.thold <- list()
  candi.pro <- list()
  candi.window.size <- list()
  tictac <- Sys.time()

  message("stage 2 of 3, evaluate individual candidate ...")

  pb <- utils::txtProgressBar(min = 0, max = n.window.size * n.pos, style = 3)
  on.exit(close(pb))
  on.exit(beepr::beep(), TRUE)

  for (i in 1:n.window.size) {
    candi.score[[i]] <- rep(0, n.pos)
    candi.thold[[i]] <- rep(0, n.pos)
    candi.pro[[i]] <- list()
    candi.window.size[[i]] <- rep(1, n.pos) * window.size[i]

    pre <- mass.pre(data, data.size, window.size = window.size[i])

    for (j in 1:n.pos) {
      dist.pro <- mass(pre$data.fft, candi[[i]][[j]], data.size, window.size[i], pre$data.mean, pre$data.sd, mean(candi[[i]][[j]]), std(candi[[i]][[j]]))
      dist.pro <- Re(sqrt(dist.pro$distance.profile))
      candi.pro[[i]][[j]] <- dist.pro
      exc.st <- max(1, candi.idx[[i]][j] - window.size[i])
      exc.ed <- min(length(dist.pro), candi.idx[[i]][j] + window.size[i])
      dist.pro[exc.st:exc.ed] <- Inf

      golden <- golden.section(dist.pro, label, pos.st, pos.ed, beta, window.size[i])

      candi.thold[[i]][j] <- golden$thold
      candi.score[[i]][j] <- golden$score

      utils::setTxtProgressBar(pb, ((i - 1) * n.pos + j))
    }
  }

  tictac <- Sys.time() - tictac
  message(sprintf("\nFinished in %.2f %s", tictac, units(tictac)))

  candi.pro.exp <- list()
  candi.exp <- list()

  for (i in 1:n.window.size) {
    candi.pro.exp[((i - 1) * n.pos + 1):(i * n.pos)] <- candi.pro[[i]]
    candi.exp[((i - 1) * n.pos + 1):(i * n.pos)] <- candi[[i]]
  }

  candi.pro <- candi.pro.exp
  candi <- candi.exp
  candi.score <- unlist(candi.score)
  candi.thold <- unlist(candi.thold)
  candi.idx <- unlist(candi.idx)
  candi.window.size <- unlist(candi.window.size)
  candi.score.sorted <- sort(signif(candi.score, 6), decreasing = TRUE, index.return = TRUE)
  candi.score <- candi.score.sorted$x
  order <- candi.score.sorted$ix

  candi.thold <- candi.thold[order]
  candi.idx <- candi.idx[order]
  candi.window.size <- candi.window.size[order]
  candi.pro <- candi.pro[order]
  candi <- candi[order]

  ## check max pattern allowed
  pat.max <- min(pat.max, floor(n.pos * 0.5))
  if (pat.max < 2) {
    return(list(score = candi.score[1], score.hist = candi.score[1], pattern = candi[[1]], thold = candi.thold[1]))
  }

  ## check combined pattern
  max.window.size <- max(window.size)
  max.pro.len <- length(data) - min(window.size) + 1
  best.pat <- rep(FALSE, n.pos * n.window.size)
  best.score <- -Inf
  exc.mask <- rep(FALSE, max.pro.len)
  score.hist <- rep(Inf, n.pos * n.window.size)
  tictac <- Sys.time()

  message("stage 3 of 3, evaluate combination of candidates ...")

  close(pb)
  pb <- utils::txtProgressBar(min = 0, max = pat.max * n.window.size * n.pos, style = 3)

  for (i in 1:pat.max) {
    pat.score <- rep(-Inf, n.pos * n.window.size)
    exc.mask.cur <- exc.mask
    exc.st <- rep(0, n.pos * n.window.size)
    exc.ed <- rep(0, n.pos * n.window.size)
    thold.cur <- list()

    for (j in 1:(n.pos * n.window.size)) {
      if (best.pat[j]) {
        next
      }

      best.pat.cur <- best.pat
      best.pat.cur[j] <- TRUE

      exc.st[j] <- max(1, candi.idx[j] - max.window.size)
      exc.ed[j] <- min(max.pro.len, candi.idx[j] + max.window.size)
      exc.mask.cur[exc.st[j]:exc.ed[j]] <- TRUE

      pro.cur <- candi.pro[best.pat.cur]
      pro.max <- -Inf
      pro.min <- Inf

      for (k in 1:length(pro.cur)) {
        pro.max <- max(max(pro.cur[[k]][!is.infinite(pro.cur[[k]])]), pro.max)
        pro.min <- min(min(pro.cur[[k]]), pro.min)
        pro.cur[[k]][exc.mask.cur] <- Inf
      }

      thold.cur[[j]] <- candi.thold[best.pat.cur]
      window.size.cur <- candi.window.size[best.pat.cur]

      iter <- 0
      score <- NULL
      while (TRUE) {
        iter <- iter + 1
        thold.old <- thold.cur[[j]]
        for (k in length(thold.cur[[j]]):1) {
          gold <- golden.section.2(
            pro.cur,
            thold.cur[[j]],
            label,
            pos.st,
            pos.ed,
            beta,
            window.size.cur[k],
            k
          )
          thold.cur[[j]] <- gold$thold
          score <- gold$score
        }

        if ((iter > 200) || (mean(thold.cur[[j]] - thold.old) < ((pro.max - pro.min) * 0.001))) {
          break
        }
      }

      pat.score[j] <- score
      exc.mask.cur[exc.st[j]:exc.ed[j]] <- FALSE

      utils::setTxtProgressBar(pb, ((i - 1) * (n.pos * n.window.size) + j))
    }

    utils::setTxtProgressBar(pb, ((i - 1) * (n.pos * n.window.size) + (n.pos * n.window.size)))

    best.candi.idx <- which.max(pat.score)

    if ((pat.score[best.candi.idx] - best.score) > 0) {
      score.hist[i] <- pat.score[best.candi.idx]
      best.score <- pat.score[best.candi.idx]
      best.pat[best.candi.idx] <- TRUE
      candi.thold[best.pat] <- thold.cur[[best.candi.idx]]
      exc.mask[exc.st[best.candi.idx]:exc.ed[best.candi.idx]] <- TRUE
    } else {
      break
    }
  }

  utils::setTxtProgressBar(pb, (pat.max * n.pos * n.window.size))

  tictac <- Sys.time() - tictac
  message(sprintf("\nFinished in %.2f %s", tictac, units(tictac)))

  score.hist <- score.hist[!is.infinite(score.hist)]

  return(list(score = best.score, score.hist = score.hist, pattern = candi[best.pat], thold = candi.thold[best.pat]))
}

#' Computes the golden section for individual candidates
#'
#' @param dist.pro the candidate distance profile
#' @param label a vector with the data bool annotation
#' @param pos.st a vector with the starting points of label
#' @param pos.ed a vector with the ending points of label
#' @param beta a number that balance the F-Score. Beta > 1 towards recall, < towards precision
#' @param window.size an integer with the sliding window size
#'
#' @return Returns the best threashold and its F-Score
#'
#' @keywords internal
#'
golden.section <- function(dist.pro, label, pos.st, pos.ed, beta, window.size) {
  golden.ratio <- (1 + sqrt(5)) / 2
  a.thold <- min(dist.pro)
  b.thold <- max(dist.pro[!is.infinite(dist.pro)])
  c.thold <- b.thold - (b.thold - a.thold) / golden.ratio
  d.thold <- a.thold + (b.thold - a.thold) / golden.ratio
  tol <- max((b.thold - a.thold) * 0.001, 0.0001)

  while (abs(c.thold - d.thold) > tol) {
    c.score <- compute.f.meas(label, pos.st, pos.ed, dist.pro, c.thold, window.size, beta)
    d.score <- compute.f.meas(label, pos.st, pos.ed, dist.pro, d.thold, window.size, beta)

    if (c.score$f.meas > d.score$f.meas) {
      b.thold <- d.thold
    } else {
      a.thold <- c.thold
    }

    c.thold <- b.thold - (b.thold - a.thold) / golden.ratio
    d.thold <- a.thold + (b.thold - a.thold) / golden.ratio
  }
  thold <- (a.thold + b.thold) * 0.5
  score <- compute.f.meas(label, pos.st, pos.ed, dist.pro, thold, window.size, beta)

  return(list(thold = thold, score = score$f.meas))
}

#' Computes the golden section for combined candidates

#' @param dist.pro the candidate distance profile
#' @param thold a number with the threshold used to calculate the F-Score
#' @param label a vector with the data bool annotation
#' @param pos.st a vector with the starting points of label
#' @param pos.ed a vector with the ending points of label
#' @param beta a number that balance the F-Score. Beta > 1 towards recall, < towards precision
#' @param window.size an integer with the sliding window size
#' @param fit.idx an integer with the index of the current threshold
#'
#' @return Returns the best threashold and its F-Score
#'
#' @keywords internal

golden.section.2 <- function(dist.pro, thold, label, pos.st, pos.ed, beta, window.size, fit.idx) {
  golden.ratio <- (1 + sqrt(5)) / 2
  a.thold <- min(dist.pro[[fit.idx]], na.rm = TRUE) ## TODO: check why NA in dist.pro
  b.thold <- max(dist.pro[[fit.idx]][!is.infinite(dist.pro[[fit.idx]])], na.rm = TRUE)
  c.thold <- b.thold - (b.thold - a.thold) / golden.ratio
  d.thold <- a.thold + (b.thold - a.thold) / golden.ratio
  tol <- max((b.thold - a.thold) * 0.001, 0.0001)

  while (abs(c.thold - d.thold) > tol) {
    c.thold.combined <- thold
    d.thold.combined <- thold
    c.thold.combined[fit.idx] <- c.thold
    d.thold.combined[fit.idx] <- d.thold

    c.score <- compute.f.meas(label, pos.st, pos.ed, dist.pro, c.thold.combined, window.size, beta)
    d.score <- compute.f.meas(label, pos.st, pos.ed, dist.pro, d.thold.combined, window.size, beta)

    if (c.score$f.meas > d.score$f.meas) {
      b.thold <- d.thold
    } else {
      a.thold <- c.thold
    }

    c.thold <- b.thold - (b.thold - a.thold) / golden.ratio
    d.thold <- a.thold + (b.thold - a.thold) / golden.ratio
  }
  thold[fit.idx] <- (a.thold + b.thold) * 0.5
  score <- compute.f.meas(label, pos.st, pos.ed, dist.pro, thold, window.size, beta)

  # beta = 2;   emphacise recall
  # beta = 0.5; emphacise precision
  return(list(thold = thold, score = score$f.meas))
}

#' Computes de F-Score
#'
#' @param label a vector with the data bool annotation
#' @param pos.st a vector with the starting points of label
#' @param pos.ed a vector with the ending points of label
#' @param dist.pro the distance profile of the data
#' @param thold a number with the threshold used to compute the prediction
#' @param window.size an integer with the sliding window size
#' @param beta a number that balance the F-Score. Beta > 1 towards recall, < towards precision
#'
#' @return Returns the F-Score, precision and recall values
#'
#' @keywords internal

compute.f.meas <- function(label, pos.st, pos.ed, dist.pro, thold, window.size, beta) {
  # generate annotation curve for each pattern
  if (is.list(dist.pro)) {
    anno.st <- list()
    n.pat <- length(dist.pro)

    for (i in 1:n.pat) {
      annor <- dist.pro[[i]] - thold[i]
      annor[annor > 0] <- 0
      annor[annor < 0] <- -1
      annor <- -annor
      anno.st[[i]] <- which(diff(c(0, annor, 0)) == 1) + 1
      anno.st[[i]] <- anno.st[[i]] - 1
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

    anno.ed <- anno.st + window.size - 1
  } else {
    anno <- dist.pro - thold
    anno[anno > 0] <- 0
    anno[anno < 0] <- -1
    anno <- -anno

    anno.st <- which(diff(c(0, anno, 0)) == 1) + 1
    anno.ed <- anno.st + window.size - 1
    anno.st <- anno.st - 1
    anno.ed <- anno.ed - 1
  }

  anno <- rep(FALSE, length(label))

  for (i in 1:length(anno.st)) {
    anno[anno.st[i]:anno.ed[i]] <- 1
  }

  is.tp <- rep(FALSE, length(anno.st))

  for (i in 1:length(anno.st)) {
    if (anno.ed[i] > length(label)) {
      anno.ed[i] <- length(label)
    }
    if (sum(label[anno.st[i]:anno.ed[i]]) > (0.8 * window.size)) {
      is.tp[i] <- TRUE
    }
  }
  tp.pre <- sum(is.tp)

  is.tp <- rep(FALSE, length(pos.st))
  for (i in 1:length(pos.st)) {
    if (sum(anno[pos.st[i]:pos.ed[i]]) > (0.8 * window.size)) {
      is.tp[i] <- TRUE
    }
  }
  tp.rec <- sum(is.tp)

  pre <- tp.pre / length(anno.st)
  rec <- tp.rec / length(pos.st)

  f.meas <- (1 + beta^2) * (pre * rec) / ((beta^2) * pre + rec)
  if (is.na(f.meas)) {
    f.meas <- 0
  }
  return(list(f.meas = f.meas, pre = pre, rec = rec))
}
