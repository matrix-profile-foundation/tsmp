#' Plot arcs between indexes of a Profile Index
#'
#' Sometimes may be useful to see where is the nearest neighbor graphically. This is the reasoning
#' behind, for example, FLUSS which uses the arc count to infer a semantic change, and SiMPle which
#' infer that arcs connect similar segments of a music. See details for a deeper explanation how to
#' use this function.
#'
#' @details
#' You have two options to use this function. First you can provide just the data, and the function
#' will try its best to retrieve the pairs for plotting. Second, you can skip the first parameters
#' and just provide the `pairs`, which is a `matrix` with two columns; the first is the starting
#' index, the second is the end index. Two colors are used to allow you to identify the direction of
#' the arc. If you use the `rpi` or `lpi` as input, you will see that these profile indexes have
#' just one direction.
#'
#' `exclusion_zone` is used to filter out small arcs that may be useless (e.g. you may be interested
#' in similarities that are far away). `edge_limit` is used to filter out spurious arcs that are
#' used connect the beginning and the end of the profile (e.g. silent audio). `threshold` is used to
#' filter indexes that have distant nearest neighbor (e.g. retrieve only the best motifs).
#'
#' @param pairs a `matrix` with 2 columns.
#' @param alpha a `numeric`. (Default is `NULL`, automatic). Alpha value for lines transparency.
#' @param quality an `int`. (Default is `30`). Number of segments to draw the arc. Bigger value,
#'   harder to render.
#' @param lwd an `int`. (Default is `15`). Line width.
#' @param col a `vector` of colors. (Default is `c("blue", "orange")`). Colors for right and left
#'   arc, respectively. Accepts one color.
#' @param main a `string`. (Default is `"Arc Plot"`). Main title.
#' @param ylab a `string`. (Default is `""`). Y label.
#' @param xlab a `string`. (Default is `"Profile Index"`). X label.
#' @param \dots further arguments to be passed to [plot()]. See [par()].
#' @param xmin an `int`. (Default is `NULL`). Set the minimum value of x axis.
#' @param xmax an `int`. (Default is `NULL`). Set the maximum value of x axis.
#'
#' @return None
#' @keywords hplot
#'
#' @export
#' @examples
#' plot_arcs(pairs = matrix(c(5, 10, 1, 10, 20, 5), ncol = 2, byrow = TRUE))
plot_arcs <- function(pairs, alpha = NULL, quality = 30, lwd = 15, col = c("blue", "orange"),
                      main = "Arc Plot", ylab = "", xlab = "Profile Index", xmin = NULL, xmax = NULL, ...) {
  if (length(pairs) == 0) {
    warning("No arc to plot.")
    return(NULL)
  }

  if (is.null(xmin)) {
    xmin <- min(pairs)
  }

  if (is.null(xmax)) {
    xmax <- max(pairs)
  }

  max_arc <- max(abs(pairs[, 2] - pairs[, 1]))
  ymax <- (max_arc / 2 + (lwd * lwd) / 8)
  z_seq <- seq(0, base::pi, length.out = quality)
  xlim <- c(xmin, xmax)
  ylim <- c(0, ymax)

  if (is.null(alpha)) {
    alpha <- min(0.5, max(10 / nrow(pairs), 0.03))
  }

  arccolr <- grDevices::adjustcolor(col, alpha.f = alpha)
  if (length(col) > 1) {
    arccoll <- grDevices::adjustcolor(col[2], alpha.f = alpha)
  } else {
    arccoll <- grDevices::adjustcolor(col, alpha.f = alpha)
  }

  # blank plot
  graphics::plot(0.5, 0.5,
    type = "n", main = main, xlab = xlab, ylab = ylab,
    xlim = xlim, ylim = ylim, yaxt = "n", ...
  )

  for (i in seq_len(nrow(pairs))) {
    if (pairs[i, 1] > pairs[i, 2]) {
      arccol <- arccoll
    } else {
      arccol <- arccolr
    }

    x1 <- min(pairs[i, 1], pairs[i, 2])
    x2 <- max(pairs[i, 1], pairs[i, 2])
    center <- (x1 - x2) / 2 + x2
    radius <- (x2 - x1) / 2
    x_seq <- center + radius * cos(z_seq)
    y_seq <- radius * sin(z_seq)
    graphics::lines(x_seq, y_seq,
      col = arccol, lwd = lwd, lty = 1, lend = 1
    )
  }

  graphics::legend(xmin, ymax,
    legend = c("Right", "Left"),
    col = grDevices::adjustcolor(col, alpha.f = 0.5), lty = 1, cex = 0.8, lwd = 5
  )
}

#' Plot a TSMP object
#'
#' @param x a Matrix Profile
#' @param data the data used to build the Matrix Profile, if not embedded to it.
#' @param type "data" or "matrix". Choose what will be plotted.
#' @param exclusion_zone if a `number` will be used instead of Matrix Profile's. (Default is `NULL`).
#' @param edge_limit if a `number` will be used instead of Matrix Profile's exclusion zone. (Default is `NULL`).
#' @param threshold the maximum value to be used to plot.
#' @param main a `string`. Main title.
#' @param xlab a `string`. X label.
#' @param ylab a `string`. Y label.
#' @param ncol an `int`. Number of columns to plot Motifs.
#' @param \dots further arguments to be passed to [plot()]. See [par()].
#'
#' @return None
#'
#' @export
#' @keywords hplot
#' @name plot
#'
#' @examples
#'
#' mp <- tsmp(mp_toy_data$data[1:200, 1], window_size = 30, verbose = 0)
#' plot(mp)
plot.ArcCount <- function(x, data, type = c("data", "matrix"), exclusion_zone = NULL, edge_limit = NULL,
                          threshold = stats::quantile(x$cac, 0.1), main = "Arcs Discover", xlab = "index",
                          ylab = "", ...) {
  def_par <- graphics::par(no.readonly = TRUE)

  if (missing(data) && !is.null(x$data)) {
    data <- x$data[[1]]
  } else {
    is.null(data) # check data presence before plotting anything
  }

  type <- match.arg(type)

  if (is.null(exclusion_zone)) {
    if (floor(x$ez * 10) < (length(x$mp) / 3)) {
      exclusion_zone <- floor(x$ez * 10)
    }
  }

  if (is.null(edge_limit)) {
    if (floor(x$ez * 10) < (length(x$mp) / 3)) {
      edge_limit <- floor(x$ez * 10)
    }
  }

  if (type == "data") {
    plot_data <- data
    data_lab <- ylab
    data_main <- "Data"
  } else {
    plot_data <- c(x$mp, rep(NA, x$w - 1))
    data_lab <- "distance"
    data_main <- "Matrix Profile"
  }

  cac <- x$cac # keep cac intact
  cac_size <- length(cac)
  profile_index <- x$pi

  if (cac_size < nrow(profile_index)) {
    warning("cac_size < profile_index")
    cac_offset <- nrow(profile_index) - cac_size
    plot_data <- as.matrix(utils::tail(plot_data, cac_size))
    profile_index <- as.matrix(utils::tail(profile_index, cac_size) - cac_offset)
  }

  xnum <- seq_len(nrow(x$mp) + min(x$w) - 1)

  offset <- attr(x, "offset")
  offset <- ifelse(is.null(offset), 0, offset)

  if (!is.null(offset)) {
    xnum <- xnum + offset
  }

  pairs <- matrix(0, cac_size, 2)
  pairs[, 1] <- seq_len(cac_size) + offset
  pairs[, 2] <- profile_index + offset

  if (threshold < min(cac)) {
    stop(paste0("`threshold` is too small for this Arc Count. Min: ", round(min(cac), 2), ", Max: ", round(max(cac), 2)))
  }

  # remove excess of arcs
  if (floor(x$w * exclusion_zone) < nrow(x$mp) / 3) {
    exclusion_zone <- floor(x$w * exclusion_zone)
  }
  if (floor(x$w * edge_limit) < length(x$mp) / 3) {
    edge_limit <- floor(x$w * edge_limit)
  }
  cac[1:edge_limit] <- Inf
  cac[(cac_size - edge_limit + 1):cac_size] <- Inf

  ind <- which(cac <= threshold)
  pairs <- pairs[ind, , drop = FALSE]

  pairdiff <- pairs[, 1] - pairs[, 2]
  ind <- which(abs(pairdiff) > exclusion_zone)
  pairs <- pairs[ind, , drop = FALSE]

  xmin <- min(xnum)
  xmax <- max(xnum)
  xlim <- c(xmin, xmax)

  graphics::layout(matrix(c(1, 2, 3), ncol = 1, byrow = TRUE))
  graphics::par(oma = c(1, 1, 3, 0), cex.lab = 1.5)
  plot_arcs(pairs, xlab = xlab, xmin = xmin, xmax = xmax, ...)
  graphics::mtext(text = main, font = 2, cex = 1.5, outer = TRUE)
  graphics::plot(xnum, c(cac, rep(NA, x$w - 1)), main = "Arc count", type = "l", xlab = xlab, ylab = "normalized count", xlim = xlim, ...)
  graphics::plot(xnum, plot_data, main = data_main, type = "l", xlab = xlab, ylab = data_lab, xlim = xlim, ...)

  graphics::par(def_par)
}

#' @export
#' @keywords hplot
#' @name plot
#'
plot.Valmod <- function(x, ylab = "distance", xlab = "index", main = "Valmod Matrix Profile", data = FALSE, ...) {
  def_par <- graphics::par(no.readonly = TRUE)
  allmatrix <- FALSE
  num_charts <- 1

  xnum <- seq_len(nrow(x$mp) + min(x$w) - 1)

  if (!is.null(attr(x, "offset"))) {
    xnum <- xnum + attr(x, "offset")
  }

  if (data) {
    num_charts <- num_charts + 1
  }

  if (num_charts > 1) {
    graphics::layout(matrix(seq_len(num_charts), ncol = 1, byrow = TRUE))
  }
  graphics::par(
    mar = c(4.1, 4.1, 2.1, 2.1),
    oma = c(1, 1, 3, 0), cex.lab = 1.5
  )

  if (data) {
    graphics::plot(xnum, x$data[[1]], type = "l", main = paste0("Data"), ylab = ylab, xlab = xlab, ...)
    graphics::mtext(text = main, font = 2, cex = 1.5, outer = TRUE)
  }

  graphics::plot(xnum, c(x$mp, rep(NA, min(x$w) - 1)), type = "l", main = paste0("Matrix Profile (w = ", min(x$w), "-", max(x$w), "; ez = ", x$ez, ")"), ylab = ylab, xlab = xlab, ...)
  graphics::mtext(text = main, font = 2, cex = 1.5, outer = TRUE)

  if (allmatrix == TRUE) {
    graphics::plot(xnum, c(x$rmp, rep(NA, min(x$w) - 1)), type = "l", main = "Right Matrix Profile", ylab = ylab, xlab = xlab, ...)
    graphics::plot(xnum, c(x$lmp, rep(NA, min(x$w) - 1)), type = "l", main = "Left Matrix Profile", ylab = ylab, xlab = xlab, ...)
  }

  graphics::par(def_par)
}

#' @export
#' @keywords hplot
#' @name plot
#'
plot.MatrixProfile <- function(x, ylab = "distance", xlab = "index", main = "Unidimensional Matrix Profile", data = FALSE, ...) {
  def_par <- graphics::par(no.readonly = TRUE)
  allmatrix <- FALSE
  num_charts <- 1

  xnum <- seq_len(nrow(x$mp) + x$w - 1)

  if (!is.null(attr(x, "offset"))) {
    xnum <- xnum + attr(x, "offset")
  }

  if (!is.null(attr(x, "join")) && !attr(x, "join")) {
    if (!any(is.null(x$rmp), is.null(x$lmp))) {
      allmatrix <- TRUE
      num_charts <- 3
    }
  }

  if (data) {
    num_charts <- num_charts + 1
  }

  if (num_charts > 1) {
    graphics::layout(matrix(seq_len(num_charts), ncol = 1, byrow = TRUE))
  }
  graphics::par(
    mar = c(4.1, 4.1, 2.1, 2.1),
    oma = c(1, 1, 3, 0), cex.lab = 1.5
  )

  if (data) {
    graphics::plot(xnum, x$data[[1]], type = "l", main = paste0("Data"), ylab = "", xlab = xlab, ...)
    graphics::mtext(text = main, font = 2, cex = 1.5, outer = TRUE)
  }

  graphics::plot(xnum, c(x$mp, rep(NA, x$w - 1)), type = "l", main = paste0("Matrix Profile (w = ", x$w, "; ez = ", x$ez, ")"), ylab = ylab, xlab = xlab, ...)
  graphics::mtext(text = main, font = 2, cex = 1.5, outer = TRUE)

  if (allmatrix == TRUE) {
    graphics::plot(xnum, c(x$rmp, rep(NA, x$w - 1)), type = "l", main = "Right Matrix Profile", ylab = ylab, xlab = xlab, ...)
    graphics::plot(xnum, c(x$lmp, rep(NA, x$w - 1)), type = "l", main = "Left Matrix Profile", ylab = ylab, xlab = xlab, ...)
  }

  graphics::par(def_par)
}

#' @export
#' @keywords hplot
#' @name plot
#'
plot.MultiMatrixProfile <- function(x, ylab = "distance", xlab = "index", main = "Multidimensional Matrix Profile", ...) {
  def_par <- graphics::par(no.readonly = TRUE)
  allmatrix <- FALSE
  num_charts <- 1
  mask <- !is.na(x$mp[1, ])
  x$mp <- x$mp[, mask, drop = FALSE]
  n_dim <- ncol(x$mp)

  xnum <- seq_len(nrow(x$mp) + x$w - 1)

  if (!is.null(attr(x, "offset"))) {
    xnum <- xnum + attr(x, "offset")
  }

  if (!is.null(attr(x, "join")) && !attr(x, "join")) {
    if (!any(is.null(x$rmp), is.null(x$lmp))) {
      allmatrix <- TRUE
      num_charts <- 3
    }
  }

  if (allmatrix == TRUE) {
    graphics::layout(matrix(seq_len(num_charts * n_dim), ncol = n_dim, byrow = TRUE))
  }

  graphics::par(
    mar = c(4.1, 4.1, 2.1, 2.1),
    oma = c(1, 1, 3, 0), cex.lab = 1.5
  )
  for (i in seq_len(n_dim)) {
    graphics::plot(xnum, c(x$mp[, i], rep(NA, min(x$w) - 1)), type = "l", main = paste0("Matrix Profile (w = ", x$w, "; ez = ", x$ez, ")"), ylab = ylab, xlab = xlab, ...)
  }
  graphics::mtext(text = main, font = 2, cex = 1.5, outer = TRUE)

  if (allmatrix == TRUE) {
    for (i in seq_len(n_dim)) {
      graphics::plot(xnum, c(x$rmp[, i], rep(NA, min(x$w) - 1)), type = "l", main = "Right Matrix Profile", ylab = ylab, xlab = xlab, ...)
    }
    for (i in seq_len(n_dim)) {
      graphics::plot(xnum, c(x$lmp[, i], rep(NA, min(x$w) - 1)), type = "l", main = "Left Matrix Profile", ylab = ylab, xlab = xlab, ...)
    }
  }

  graphics::par(def_par)
}

#' @export
#' @keywords hplot
#' @name plot
#'
plot.SimpleMatrixProfile <- function(x, ylab = "distance", xlab = "index", main = "SiMPle Matrix Profile", data = FALSE, ...) {
  def_par <- graphics::par(no.readonly = TRUE)
  num_charts <- 1

  xnum <- seq_len(nrow(x$mp) + x$w - 1)

  if (!is.null(attr(x, "offset"))) {
    xnum <- xnum + attr(x, "offset")
  }

  if (data) {
    num_charts <- num_charts + length(x$data)
  }

  if (num_charts > 1) {
    graphics::layout(matrix(seq_len(num_charts), ncol = 1, byrow = TRUE))
  }

  graphics::par(
    mar = c(4.1, 4.1, 2.1, 2.1),
    oma = c(1, 1, 3, 0), cex.lab = 1.5
  )

  if (data) {
    n_dim <- ncol(x$data[[1]])

    graphics::plot(xnum, x$data[[1]][, 1], type = "l", main = paste0("Data"), ylab = "", xlab = xlab, ...)
    graphics::mtext(text = main, font = 2, cex = 1.5, outer = TRUE)

    if (n_dim > 1) {
      for (i in 2:n_dim) {
        graphics::lines(xnum, x$data[[1]][, i], main = paste0("Data"), ylab = "", xlab = xlab, col = i, ...)
      }
    }

    if (length(x$data) > 1) {
      n_dim <- ncol(x$data[[2]])

      graphics::plot(xnum, x$data[[2]][, 1], type = "l", main = paste0("Query"), ylab = "", xlab = xlab, ...)
      graphics::mtext(text = main, font = 2, cex = 1.5, outer = TRUE)

      if (n_dim > 1) {
        for (i in 2:n_dim) {
          graphics::lines(xnum, x$data[[2]][, i], main = paste0("Data"), ylab = "", xlab = xlab, col = i, ...)
        }
      }
    }
  }

  graphics::plot(xnum, c(x$mp, rep(NA, min(x$w) - 1)), type = "l", main = paste0("Matrix Profile (w = ", x$w, "; ez = ", x$ez, ")"), ylab = ylab, xlab = xlab, ...)

  graphics::mtext(text = main, font = 2, cex = 1.5, outer = TRUE)
  graphics::par(def_par)
}

#' @export
#' @keywords hplot
#' @name plot
#'
plot.Fluss <- function(x, data, type = c("data", "matrix"),
                       main = "Fast Low-cost Unipotent Semantic Segmentation", xlab = "index",
                       ylab = "", ...) {
  def_par <- graphics::par(no.readonly = TRUE)

  if (missing(data) && !is.null(x$data)) {
    data <- x$data[[1]]
  } else {
    is.null(data) # check data presence before plotting anything
  }

  type <- match.arg(type)

  if (type == "data") {
    plot_data <- data
    data_lab <- ylab
    data_main <- "Data"
  } else {
    plot_data <- c(x$mp, rep(NA, min(x$w) - 1))
    data_lab <- "distance"
    data_main <- "Matrix Profile"
  }

  fluss_idx <- sort(x$fluss)

  fluss_size <- length(fluss_idx) + 1
  pairs <- matrix(0, fluss_size, 2)

  offset <- attr(x, "offset")
  offset <- ifelse(is.null(offset), 0, offset)

  for (i in seq_len(fluss_size)) {
    if (i == 1) {
      pairs[i, 1] <- offset
    } else {
      pairs[i, 1] <- fluss_idx[i - 1] + offset
    }

    if (i == fluss_size) {
      pairs[i, 2] <- nrow(x$mp) + offset
    } else {
      pairs[i, 2] <- fluss_idx[i] + offset
    }
  }

  xnum <- seq_len(nrow(x$mp) + min(x$w) - 1) + offset
  xmin <- min(xnum)
  xmax <- max(xnum)
  xlim <- c(xmin, xmax)

  graphics::layout(matrix(c(1, 2, 3), ncol = 1, byrow = TRUE))
  graphics::par(oma = c(1, 1, 3, 0), cex.lab = 1.5)
  plot_arcs(pairs, xlab = xlab, xmin = xmin, xmax = xmax, ...)
  graphics::mtext(text = main, font = 2, cex = 1.5, outer = TRUE)
  graphics::plot(xnum, plot_data, main = data_main, type = "l", xlab = xlab, ylab = data_lab, xlim = xlim, ...)
  graphics::plot(xnum, c(x$cac, rep(NA, min(x$w) - 1)), main = "Arc count", type = "l", xlab = xlab, ylab = "normalized count", xlim = xlim, ylim = c(0, 1), ...)

  graphics::par(def_par)
}

#' @export
#' @keywords hplot
#' @name plot
#'
plot.Floss <- function(x, data, type = c("data", "matrix"),
                       main = "Fast Low-cost Online Semantic Segmentation", xlab = "index",
                       ylab = "", ...) {
  def_par <- graphics::par(no.readonly = TRUE)

  if (missing(data) && !is.null(x$data)) {
    data <- x$data[[1]]
  } else {
    is.null(data) # check data presence before plotting anything
  }

  type <- match.arg(type)

  if (type == "data") {
    plot_data <- data
    data_lab <- ylab
    data_main <- "Data"
  } else {
    plot_data <- c(x$mp, rep(NA, min(x$w) - 1))
    data_lab <- "distance"
    data_main <- "Matrix Profile"
  }

  offset <- attr(x, "offset")
  offset <- ifelse(is.null(offset), 0, offset)

  floss_idx <- sort(x$floss)
  floss_idx <- floss_idx - offset

  floss_size <- length(floss_idx) + 1
  pairs <- matrix(0, floss_size, 2)

  for (i in seq_len(floss_size)) {
    if (i == 1) {
      pairs[i, 1] <- 0 # offset
    } else {
      pairs[i, 1] <- floss_idx[i - 1] + offset
    }

    if (i == floss_size) {
      pairs[i, 2] <- nrow(x$mp) + offset
    } else {
      pairs[i, 2] <- floss_idx[i] + offset
    }
  }

  xnum <- seq_len(nrow(x$mp) + min(x$w) - 1) + offset
  xmin <- min(xnum)
  xmax <- max(xnum)
  xlim <- c(xmin, xmax)

  cac_fin_len <- length(x$cac_final)
  mp_len <- nrow(x$mp)
  new_data <- attr(x, "new_data")

  if (cac_fin_len == floor((mp_len * vars()$kmode + new_data / 2))) {
    cac <- x$cac_final
  } else {
    cac <- utils::tail(x$cac_final, -offset)
  }
  cac_size <- length(cac)
  cac <- c(cac, rep(NA, nrow(x$mp) + min(x$w) - 1 - cac_size))

  graphics::layout(matrix(c(1, 2, 3), ncol = 1, byrow = TRUE))
  graphics::par(oma = c(1, 1, 3, 0), cex.lab = 1.5)
  plot_arcs(pairs, xlab = xlab, xmin = xmin, xmax = xmax, ...)
  graphics::mtext(text = main, font = 2, cex = 1.5, outer = TRUE)
  graphics::plot(xnum, plot_data, main = data_main, type = "l", xlab = xlab, ylab = data_lab, xlim = xlim, ...)
  graphics::plot(xnum, cac, main = "Arc count", type = "l", xlab = xlab, ylab = "normalized count", xlim = xlim, ylim = c(0, 1), ...)

  graphics::par(def_par)
}

#' @export
#' @keywords hplot
#' @name plot
#'
plot.Chain <- function(x, data, type = c("data", "matrix"), main = "Chain Discover", xlab = "index", ylab = "", ...) {
  def_par <- graphics::par(no.readonly = TRUE)

  if (missing(data) && !is.null(x$data)) {
    data <- x$data[[1]]
  } else {
    is.null(data) # check data presence before plotting anything
  }

  type <- match.arg(type)

  if (type == "data") {
    plot_data <- data
    plot_subtitle <- "Data"
  } else {
    plot_data <- c(x$mp, rep(NA, min(x$w) - 1))
    ylab <- "distance"
    plot_subtitle <- paste0("Matrix Profile (w = ", x$w, "; ez = ", x$ez, ")")
  }

  chain_size <- length(x$chain$best)
  pairs <- matrix(0, chain_size - 1, 2)

  offset <- attr(x, "offset")
  offset <- ifelse(is.null(offset), 0, offset)

  for (i in seq_len(chain_size - 1)) {
    pairs[i, 1] <- x$chain$best[i] + offset
    pairs[i, 2] <- x$chain$best[i + 1] + offset
  }

  xnum <- seq_len(nrow(x$mp) + min(x$w) - 1) + offset
  xmin <- min(xnum)
  xmax <- max(xnum)
  xlim <- c(xmin, xmax)

  # plot matrix profile
  graphics::layout(matrix(c(1, 2, 3), ncol = 1, byrow = TRUE))
  graphics::par(oma = c(1, 1, 3, 0), cex.lab = 1.5)
  plot_arcs(pairs, xlab = xlab, xmin = xmin, xmax = xmax, ...)
  graphics::mtext(text = main, font = 2, cex = 1.5, outer = TRUE)
  graphics::plot(xnum, plot_data,
    type = "l", main = plot_subtitle,
    xlim = xlim, xlab = xlab, ylab = ylab, ...
  )
  graphics::abline(v = x$chain$best + offset, col = 1:chain_size, lwd = 2)

  # blank plot
  motif <- znorm(data[x$chain$best[1]:min((x$chain$best[1] + x$w - 1), nrow(data))])
  graphics::plot(motif,
    type = "l", main = "Motifs", xlab = "length", ylab = "normalized data",
    xlim = c(0, length(motif)), ylim = c(min(motif) - chain_size / 2, max(motif)), ...
  )

  for (i in 2:chain_size) {
    motif <- znorm(data[x$chain$best[i]:min((x$chain$best[i] + x$w - 1), nrow(data))])

    graphics::lines(motif - i / 2, col = i, ...)
  }

  graphics::par(def_par)
}

#' @export
#' @keywords hplot
#' @name plot
#'
plot.Discord <- function(x, data, type = c("data", "matrix"), ncol = 3, main = "Discord Discover", xlab = "index", ylab = "", ...) {
  def_par <- graphics::par(no.readonly = TRUE)

  if (missing(data) && !is.null(x$data)) {
    data <- x$data[[1]]
  } else {
    is.null(data) # check data presence before plotting anything
  }

  if ("PMP" %in% class(x)) {
    x$mp <- as.matrix(x$pmp[[1]])
    x$pi <- as.matrix(x$pmpi[[1]])
  }

  type <- match.arg(type)

  if (type == "data") {
    plot_data <- data
    plot_subtitle <- "Data"
  } else {
    plot_data <- x$mp
    ylab <- "distance"
    plot_subtitle <- paste0("Matrix Profile (w = ", x$w, "; ez = ", x$ez, ")")
  }

  discords <- x$discord$discord_idx
  n_discords <- length(x$discord$discord_idx)
  neighbors <- x$discord$discord_neighbor
  matrix_profile_size <- nrow(x$mp)

  # layout: matrix profile on top, discords below.
  graphics::layout(matrix(
    c(rep(1, ncol), (seq_len(ceiling(n_discords / ncol) * ncol) + 1)),
    ceiling(n_discords / ncol) + 1,
    ncol,
    byrow = TRUE
  ))
  # plot matrix profile
  xnum <- seq_len(nrow(x$mp) + min(x$w) - 1)

  offset <- attr(x, "offset")
  offset <- ifelse(is.null(offset), 0, offset)

  if (!is.null(offset)) {
    xnum <- xnum + offset
  }

  graphics::par(oma = c(1, 1, 3, 0), cex.lab = 1.5)
  graphics::plot(xnum, plot_data, type = "l", main = plot_subtitle, xlab = xlab, ylab = ylab)
  graphics::mtext(text = main, font = 2, cex = 1.5, outer = TRUE)
  graphics::abline(v = unlist(discords) + offset, col = seq_len(n_discords), lwd = 3)
  graphics::abline(v = unlist(neighbors) + offset, col = rep(seq_len(n_discords), sapply(neighbors, length)), lwd = 1, lty = 2)
  # plot discords
  for (i in 1:n_discords) {
    discord1 <- znorm(data[discords[[i]]:min((discords[[i]] + x$w - 1), matrix_profile_size)])

    # blank plot
    graphics::plot(0.5, 0.5,
      type = "n", main = paste("Discord", i), xlab = "length", ylab = "normalized data",
      xlim = c(0, length(discord1)), ylim = c(min(discord1), max(discord1))
    )

    for (j in seq_len(length(neighbors[[i]]))) {
      neigh <- znorm(data[neighbors[[i]][j]:min((neighbors[[i]][j] + x$w - 1), matrix_profile_size)])
      graphics::lines(neigh, col = "gray70", lty = 2)
    }

    graphics::lines(discord1, col = i, lwd = 2)
  }

  graphics::par(def_par)
}


#' @export
#' @keywords hplot
#' @name plot
#'

plot.Snippet <- function(x, data, ncol = 3, main = "Snippet Finder", xlab = "index", ylab = "", ...) {
  def_par <- graphics::par(no.readonly = TRUE)

  snippets <- x$snippet_idx
  n_snippets <- length(x$snippet_idx)

  if (n_snippets == 0) {
    graphics::par(def_par)
    stop("No Snippets found to plot.")
  }

  if (missing(data) && !is.null(x$data)) {
    data <- x$data[[1]]
  } else {
    is.null(data) # check data presence before plotting anything
  }

  plot_data <- data
  plot_subtitle <- "Data"

  # layout: matrix profile on top, motifs below.
  graphics::layout(matrix(
    c(rep(1, ncol), rep(2, ncol), (seq_len(ceiling(n_snippets / ncol) * ncol) + 2)),
    ceiling(n_snippets / ncol) + 2,
    ncol,
    byrow = TRUE
  ))

  # plot data
  xnum <- seq_len(nrow(data))

  offset <- attr(x, "offset")
  offset <- ifelse(is.null(offset), 0, offset)

  if (!is.null(offset)) {
    xnum <- xnum + offset
  }

  graphics::par(oma = c(1, 1, 3, 0), cex.lab = 1.5)
  graphics::plot(xnum, plot_data, type = "l", main = plot_subtitle, xlab = xlab, ylab = ylab)
  graphics::mtext(text = main, font = 2, cex = 1.5, outer = TRUE)

  graphics::plot(xnum, rep(1, length(xnum)),
    ylab = "", xlab = "Index",
    type = "p", main = "Horizontal regime bar", pch = 15, cex = 0.5,
    col = x$regime + 1
  )

  for (i in 1:n_snippets) {
    snip <- znorm(data[snippets[i]:min((snippets[i] + x$snippet_size - 1), nrow(data))])

    # blank plot
    graphics::plot(0.5, 0.5,
      type = "n", main = paste("Snippet", i), xlab = "length", ylab = "normalized data",
      xlim = c(0, length(snip)), ylim = c(min(snip), max(snip))
    )

    graphics::lines(snip, col = i + 1, lwd = 2)
  }

  # obj <- list(snippet_idx = snippetidx, snippet_frac = fraction, snippet_size = s_size, regime = horizontal, data = list(data))


  graphics::par(def_par)
}


#' @export
#' @keywords hplot
#' @name plot
#'
plot.Motif <- function(x, data, type = c("data", "matrix"), ncol = 3, main = "MOTIF Discover", xlab = "index", ylab = "", ...) {
  def_par <- graphics::par(no.readonly = TRUE)

  motifs <- x$motif$motif_idx
  n_motifs <- length(x$motif$motif_idx)

  if (n_motifs == 0) {
    graphics::par(def_par)
    stop("No Motifs found to plot.")
  }

  if ("PMP" %in% class(x)) {
    x$mp <- as.matrix(x$pmp[[1]])
    x$pi <- as.matrix(x$pmpi[[1]])
  }

  if ("Valmod" %in% class(x)) {
    valmod <- TRUE

    if (main == "MOTIF Discover") {
      main <- paste("Valmod", main)
    }
  } else {
    valmod <- FALSE
  }

  if (missing(data) && !is.null(x$data)) {
    data <- x$data[[1]]
  } else {
    is.null(data) # check data presence before plotting anything
  }

  type <- match.arg(type)

  if (type == "data") {
    plot_data <- data
    plot_subtitle <- "Data"
  } else {
    plot_data <- c(x$mp, rep(NA, min(x$w) - 1))
    ylab <- "distance"
    if (valmod) {
      plot_subtitle <- paste0("Matrix Profile (w = ", min(x$w), "-", max(x$w), "; ez = ", x$ez, ")")
    } else {
      plot_subtitle <- paste0("Matrix Profile (w = ", min(x$w), "; ez = ", x$ez, ")")
    }
  }

  neighbors <- x$motif$motif_neighbor
  windows <- unlist(x$motif$motif_window)

  # layout: matrix profile on top, motifs below.
  graphics::layout(matrix(
    c(rep(1, ncol), (seq_len(ceiling(n_motifs / ncol) * ncol) + 1)),
    ceiling(n_motifs / ncol) + 1,
    ncol,
    byrow = TRUE
  ))
  # plot matrix profile
  xnum <- seq_len(nrow(x$mp) + min(x$w) - 1)

  offset <- attr(x, "offset")
  offset <- ifelse(is.null(offset), 0, offset)

  if (!is.null(offset)) {
    xnum <- xnum + offset
  }
  graphics::par(oma = c(1, 1, 3, 0), cex.lab = 1.5)
  graphics::plot(xnum, plot_data, type = "l", main = plot_subtitle, xlab = xlab, ylab = ylab)
  graphics::mtext(text = main, font = 2, cex = 1.5, outer = TRUE)
  graphics::abline(v = unlist(motifs) + offset, col = rep(1:n_motifs, each = 2), lwd = c(3, 1))
  graphics::abline(v = unlist(neighbors) + offset, col = rep(1:n_motifs, sapply(neighbors, length)), lwd = 1, lty = 2)

  # plot motifs
  if (valmod) {
    for (i in 1:n_motifs) {
      motif1 <- znorm(data[motifs[[i]][1]:min((motifs[[i]][1] + windows[i] - 1), nrow(data))])
      motif2 <- znorm(data[motifs[[i]][2]:min((motifs[[i]][2] + windows[i] - 1), nrow(data))])

      # blank plot
      graphics::plot(0.5, 0.5,
        type = "n", main = paste0("Motif ", i, " (w = ", windows[i], ")"), xlab = "length", ylab = "normalized data",
        xlim = c(0, length(motif1)), ylim = c(min(motif1), max(motif1))
      )

      for (j in seq_len(length(neighbors[[i]]))) {
        neigh <- znorm(data[neighbors[[i]][j]:min((neighbors[[i]][j] + windows[i] - 1), nrow(data))])
        graphics::lines(neigh, col = "gray70", lty = 2)
      }

      graphics::lines(motif2, col = "black")
      graphics::lines(motif1, col = i, lwd = 2)
    }
  } else {
    for (i in 1:n_motifs) {
      motif1 <- znorm(data[motifs[[i]][1]:min((motifs[[i]][1] + x$w - 1), nrow(data))])
      motif2 <- znorm(data[motifs[[i]][2]:min((motifs[[i]][2] + x$w - 1), nrow(data))])

      # blank plot
      graphics::plot(0.5, 0.5,
        type = "n", main = paste("Motif", i), xlab = "length", ylab = "normalized data",
        xlim = c(0, length(motif1)), ylim = c(min(motif1), max(motif1))
      )

      for (j in seq_len(length(neighbors[[i]]))) {
        neigh <- znorm(data[neighbors[[i]][j]:min((neighbors[[i]][j] + x$w - 1), nrow(data))])
        graphics::lines(neigh, col = "gray70", lty = 2)
      }

      graphics::lines(motif2, col = "black")
      graphics::lines(motif1, col = i, lwd = 2)
    }
  }

  graphics::par(def_par)
}

#' @export
#' @keywords hplot
#' @name plot
#'
plot.MultiMotif <- function(x, data, type = c("data", "matrix"), ncol = 3, main = "Multidimensional MOTIF Discover", xlab = "index", ylab = "", ...) {
  def_par <- graphics::par(no.readonly = TRUE)

  motifs <- x$motif$motif_idx
  n_motifs <- length(x$motif$motif_idx)

  if (n_motifs == 0) {
    graphics::par(def_par)
    stop("No Motifs found to plot.")
  }

  if (missing(data) && !is.null(x$data)) {
    data <- x$data[[1]]
  } else {
    is.null(data) # check data presence before plotting anything
  }

  type <- match.arg(type)

  if (type == "data") {
    plot_data <- data
    plot_subtitle <- "Data"
  } else {
    plot_data <- apply(x$mp, 2, function(y) c(y, rep(NA, min(x$w) - 1)))
    ylab <- "distance"
  }

  n_dim <- x$n_dim
  motifs_dim <- x$motif$motif_dim

  dim_idx <- list()
  for (i in seq_len(n_dim)) {
    mot <- vector(mode = "numeric")
    for (j in seq_len(n_motifs)) {
      if (i %in% motifs_dim[[j]]) {
        mot <- c(mot, j)
        dim_idx[[i]] <- mot
      }
    }
  }

  # layout: matrix profile on top, motifs below.
  graphics::layout(matrix(
    c(rep(seq_len(n_dim), each = ncol), (seq_len(ceiling(n_motifs / ncol) * ncol) + n_dim)),
    # ceiling(n_motifs / ncol) + 1,
    ncol = ncol,
    byrow = TRUE
  ))
  # plot matrix profile
  xnum <- seq_len(nrow(x$mp) + min(x$w) - 1)

  offset <- attr(x, "offset")
  offset <- ifelse(is.null(offset), 0, offset)

  if (!is.null(offset)) {
    xnum <- xnum + offset
  }

  graphics::par(
    mar = c(4.1, 4.1, 2.1, 2.1),
    oma = c(1, 1, 3, 0), cex.lab = 1.5
  )
  for (i in seq_len(length(dim_idx))) {
    if (type == "matrix") {
      plot_subtitle <- paste0("Matrix Profile ", i, "d (w = ", x$w, "; ez = ", x$ez, ")")
    }

    graphics::plot(xnum, plot_data[, i],
      type = "l",
      main = plot_subtitle,
      xlab = xlab, ylab = ylab
    )

    midx <- dim_idx[[i]]
    if (!is.null(midx)) {
      graphics::abline(v = unlist(motifs[midx]) + offset, col = rep(midx, each = 2), lwd = c(2, 1))
    }
  }

  graphics::mtext(text = main, font = 2, cex = 1.5, outer = TRUE)

  # plot motifs
  for (i in 1:n_motifs) {
    dim_len <- length(motifs_dim[[i]])

    motif1 <- list()
    motif2 <- list()

    for (j in seq_len(dim_len)) {
      motif1[[j]] <- znorm(data[motifs[[i]][1]:min((motifs[[i]][1] + x$w - 1), nrow(data)), motifs_dim[[i]][j]])
      motif2[[j]] <- znorm(data[motifs[[i]][2]:min((motifs[[i]][2] + x$w - 1), nrow(data)), motifs_dim[[i]][j]])
    }


    # blank plot
    graphics::plot(0.5, 0.5,
      type = "n", main = paste("Motif", i), xlab = "length", ylab = "normalized data",
      xlim = c(0, length(motif1[[1]])), ylim = c(min(unlist(motif1), unlist(motif2)), max(unlist(motif1), unlist(motif2)))
    )

    if (length(motif2) > 1) {
      for (j in (seq_len(dim_len - 1) + 1)) {
        graphics::lines(motif2[[j]], col = i, lwd = 1)
        graphics::lines(motif1[[j]], col = i, lwd = 2)
      }
    }

    graphics::lines(motif2[[1]], col = i, lwd = 1)
    graphics::lines(motif1[[1]], col = i, lwd = 2)
  }

  graphics::par(def_par)
}


#' @export
#' @keywords hplot
#' @name plot
#'

plot.Salient <- function(x, data, main = "Salient Subsections", xlab = "index", ylab = "", ...) {
  def_par <- graphics::par(no.readonly = TRUE)

  if (missing(data) && !is.null(x$data)) {
    data <- x$data[[1]]
  } else {
    is.null(data) # check data presence before plotting anything
  }

  plot_data <- data
  plot_subtitle <- "Data"
  y_min <- min(data)
  y_max <- max(data)


  offset <- attr(x, "offset")
  offset <- ifelse(is.null(offset), 0, offset)

  mds <- salient_mds(x, data)
  idxs <- sort(x$salient$indexes[, 1]) + offset

  # layout: matrix profile on top, motifs below.
  graphics::layout(matrix(c(1, 1, 1, 0, 2, 0), ncol = 3, byrow = TRUE))
  # plot matrix profile
  xnum <- seq_len(nrow(x$mp) + min(x$w) - 1)

  if (!is.null(offset)) {
    xnum <- xnum + offset
  }
  graphics::par(oma = c(1, 1, 3, 0), cex.lab = 1.5)
  graphics::plot(xnum, plot_data, type = "l", main = plot_subtitle, xlab = xlab, ylab = ylab)
  graphics::mtext(text = main, font = 2, cex = 1.5, outer = TRUE)

  graphics::rect(idxs, y_min,
    xright = idxs + x$w, y_max, border = NA,
    col = grDevices::adjustcolor("blue", alpha.f = 0.1)
  )

  graphics::plot(mds, main = "MDS")

  graphics::par(def_par)
}

skimp_plot_set_canvas <- function(..., pmp_obj = NULL) {
  if (!is.null(pmp_obj)) {
    xmin <- 1
    ymin <- min(pmp_obj$w)
    ymax <- max(pmp_obj$w) + floor((max(pmp_obj$w) - ymin) / 24) # arbitrary
    mp_min <- length(pmp_obj$pmp[[as.character(ymin)]])
    xmax <- mp_min + ymin - 1
  } else {
    pars <- list(...)
    xmin <- pars$xmin
    xmax <- pars$xmax
    ymin <- pars$ymin
    ymax <- pars$ymax

    graphics::plot(c(xmin, xmax), c(ymin, ymax),
      main = "Pan Matrix Profile", xlab = "", ylab = "window",
      type = "n", xaxt = "n", yaxt = "n", xlim = c(xmin, xmax)
    )

    graphics::axis(side = 1, at = floor(c(xmin, seq(xmin, xmax, length.out = 10), xmax))) # X
    graphics::axis(side = 2, at = floor(c(ymin, seq(ymin, ymax, length.out = 10), ymax))) # Y
  }

  invisible()
}

skimp_plot_add_layer <- function(layer, window, window_set = NULL, func = NULL) {
  coords <- graphics::par("usr")
  xmin <- 1
  ymin <- window
  data_size <- length(layer) + window - 1 # theoretical data size

  # assert
  if (data_size != (coords[[2]] + coords[[1]] - xmin)) {
    stop("data_size calc is wrong")
  }

  if (is.null(window_set)) {
    # if this is the first layer
    w_min <- window
  } else {
    # else, find the position to plot
    w_min <- min(window_set)
    ymax <- coords[[4]] + coords[[3]] - w_min
    # assert
    # if (ymax != (max(window_set) + floor((max(window_set) - min(window_set)) / 24))) {
    #   print(list(
    #     ymax = ymax,
    #     max_window = max(window_set),
    #     min_window = min(window_set),
    #     result = (max(window_set) + floor((max(window_set) - min(window_set)) / 24))
    #   ))
    #   stop("ymax calc is wrong.")
    #   print(str(ymax = ymax, window_set = window_set, max_window = max(window_set), min_window = min(window_set)))
    # }

    upper_windows <- window_set[window_set > window]

    if (length(upper_windows) > 0) {
      ytop <- min(upper_windows)
    } else {
      ytop <- ymax
    }
  }

  layer <- c(layer, rep(0, window - 1))
  #  layer <- ed_corr(layer, window)
  xmax <- length(layer)

  # layer <- normalize(layer, 0, 1)

  if (is.function(func)) {
    layer <- func(layer)
  }

  layer[layer > 1] <- 1
  layer[layer < 0] <- 0

  # print(list(
  #   coords = coords, xmin = xmin, xmax = xmax,
  #   ymin = ymin, ymax = ymax, ytop = ytop, w_min = w_min
  # ))
  #
  message("layer: ", ymin, "-", ytop)

  graphics::image(matrix(layer, nrow = 1),
    xlim = c(xmin, xmax), ylim = c(ymin, ytop)
  )

  # graphics::rasterImage(matrix(layer, nrow = 1),
  #   xleft = xmin, xright = xmax,
  #   ybottom = ymin, ytop = ytop,
  #   interpolate = FALSE
  # )

  Sys.sleep(1) # needed for plot update

  invisible()
}

# image() ?

skimp_plot_add_raster <- function(layer, window, window_set = NULL, func = NULL) {
  coords <- graphics::par("usr")
  xmin <- 1
  ymin <- window
  data_size <- length(layer) + window - 1 # theoretical data size

  # assert
  if (data_size != (coords[[2]] + coords[[1]] - xmin)) {
    stop("data_size calc is wrong")
  }

  if (is.null(window_set)) {
    # if this is the first layer
    w_min <- window
  } else {
    # else, find the position to plot
    w_min <- min(window_set)
    ymax <- coords[[4]] + coords[[3]] - w_min
    # assert
    # if (ymax != (max(window_set) + floor((max(window_set) - min(window_set)) / 24))) {
    #   print(list(
    #     ymax = ymax,
    #     max_window = max(window_set),
    #     min_window = min(window_set),
    #     result = (max(window_set) + floor((max(window_set) - min(window_set)) / 24))
    #   ))
    #   stop("ymax calc is wrong.")
    #   print(str(ymax = ymax, window_set = window_set, max_window = max(window_set), min_window = min(window_set)))
    # }

    upper_windows <- window_set[window_set > window]

    if (length(upper_windows) > 0) {
      ytop <- min(upper_windows)
    } else {
      ytop <- ymax
    }
  }

  layer <- c(layer, rep(0, window - 1))
  #  layer <- ed_corr(layer, window)
  xmax <- length(layer)

  # layer <- normalize(layer, 0, 1)

  if (is.function(func)) {
    layer <- func(layer)
  }

  layer[layer > 1] <- 1
  layer[layer < 0] <- 0

  # print(list(
  #   coords = coords, xmin = xmin, xmax = xmax,
  #   ymin = ymin, ymax = ymax, ytop = ytop, w_min = w_min
  # ))
  #
  message("layer: ", ymin, "-", ytop)

  # graphics::rasterImage(matrix(layer, nrow = 1),
  #   xleft = xmin, xright = xmax,
  #   ybottom = ymin, ytop = ytop,
  #   interpolate = FALSE
  # )

  ras <- raster::raster(matrix(layer, nrow = 1),
    xmn = xmin, xmx = xmax,
    ymn = ymin, ymx = ytop
  )

  graphics::plot(raster::brick(ras), add = TRUE)
  Sys.sleep(1) # needed for plot update

  invisible()
}

#' @export
#' @keywords hplot
#' @name plot
#'
plot.PMP <- function(x, ylab = "distance", xlab = "index", main = "Unidimensional Matrix Profile", data = FALSE, ...) {
  def_par <- graphics::par(no.readonly = TRUE)
  # prepare plot using the values in `windows` vector.
  min_window <- min(x$w)
  max_window <- max(x$w)
  max_len <- length(x$pmp[[as.character(min_window)]])
  data_size <- max_len + min_window - 1

  if (!(max_len > 0)) {
    stop("matrix profile with window size ", min_window, " is not in the object. Cannot go further.")
  }

  # sort pmp
  idxs <- as.numeric(names(x$pmp))
  idxs <- sort(idxs, index.return = T)$ix
  all_profiles <- x$pmp[idxs]

  skimp_plot_set_canvas(
    ymin = min_window,
    ymax = max_window + floor((max_window - min_window) / 24), # arbitrary
    xmin = 1,
    xmax = data_size
  )
  Sys.sleep(1) # needed for plot update

  # now start to print all layers
  for (i in seq_along(all_profiles)) {
    if (!is.null(all_profiles[i])) {
      # layer <- tsmp:::normalize(all_profiles[[i]])
      layer <- all_profiles[[i]]
      layer[layer > 1] <- 1
      curr_w <- as.numeric(names(all_profiles[i]))
      next_w <- as.numeric(names(all_profiles[i + 1]))
      next_w <- ifelse(is.na(next_w), curr_w, next_w)
      graphics::rasterImage(matrix(layer, nrow = 1),
        xleft = 1, xright = length(all_profiles[[i]]),
        ybottom = curr_w, ytop = next_w + 5
      )
    }
  }

  graphics::par(def_par)
}

#' @keywords internal
#' @noRd
#'
plot_skimp <- function(pmp, func = NULL) {
  def_par <- graphics::par(no.readonly = TRUE)
  # prepare plot using the values in `windows` vector.
  min_window <- min(pmp$w)
  max_window <- max(pmp$w)
  mp_len <- length(pmp$pmp[[as.character(min_window)]])

  if (!(mp_len > 0)) {
    stop("matrix profile with window size ", min_window, " is not in the object. Cannot go further.")
  }

  data_size <- mp_len + min_window - 1
  sizes <- sort(pmp$w, index.return = TRUE)
  window_sizes <- sizes$x
  # window_idxs <- sizes$ix

  skimp_plot_set_canvas(
    ymin = min_window,
    ymax = max_window + floor((max_window - min_window) / 24), # arbitrary
    xmin = 1,
    xmax = data_size
  )
  Sys.sleep(1) # needed for plot update

  # now start to print all layers



  for (i in window_sizes) {
    layer <- pmp$pmp[[as.character(i)]]

    if (!is.null(layer)) {
      skimp_plot_add_layer(layer, i, window_sizes, func)
    }
  }

  graphics::par(def_par)
}
