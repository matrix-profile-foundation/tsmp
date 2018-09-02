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
#' `exclusion.zone` is used to filter out small arcs that may be useless (e.g. you may be interested
#' in similarities that are far away). `edge.limit` is used to filter out spurious arcs that are
#' used connect the beginning and the end of the profile (e.g. silent audio). `threshold` is used to
#' filter indexes that have distant nearest neighbor (e.g. retrieve only the best motifs).
#'
#' @param mp a Matrix Profile, or a Corrected Arc Count profile from [fluss.cac()].
#' @param pi a Profile Index.
#' @param window.size an `int`. Size of the sliding window.
#' @param exclusion.zone an `int`. (Default is `5`). Exclusion zone for small distances. See
#'   details.
#' @param edge.limit an `int`. (Default is `5`). Exclusion zone for the edges of the profile. See
#'   details.
#' @param threshold a `numeric`. (Default is 10% of the lowest values). Threshold for retrieving
#'   pairs. See details.
#' @param pairs a `matrix` with 2 columns. (Default is `NULL`). Instead of `mp`, `pi`,
#'   `window.size`, `edge.limit`, `threshold`, you can give your custom list of pairs.
#' @param alpha a `numeric`. (Default is `NULL`, automatic). Alpha value for lines transparency.
#' @param quality an `int`. (Default is `30`). Number of segments to draw the arc. Bigger value,
#'   harder to render.
#' @param lwd an `int`. (Default is `15`). Line width.
#' @param col a `vector` of colors. (Default is `c("blue", "orange")`). Colors for right and left
#'   arc, respectively. Accepts one color.
#' @param main a `string`. (Default is `"Arc Plot"`). Main title.
#' @param ylab a `string`. (Default is `""`). Y label
#' @param xlab a `string`. (Default is `"Profile Index"`). X label.
#' @param ... further arguments to be passed to [plot()]. See [par()].
#'
#' @return Quietly returns the computed `pairs`, or those you gave as input.
#' @export
#'
#' @examples
#' arcplot(pairs = matrix(c(5, 10, 1, 10, 20, 5), ncol = 2, byrow = TRUE))
arcplot <- function(mp, pi, window.size, exclusion.zone = 5, edge.limit = 5, threshold = quantile(mp, 0.1), pairs = NULL, alpha = NULL, quality = 30,
                    lwd = 15, col = c("blue", "orange"), main = "Arc Plot", ylab = "", xlab = "Profile Index",
                    ...) {
  if (is.null(pairs)) {
    data <- mp
    data.size <- nrow(data)
    pairs <- matrix(0, nrow(mp), 2)
    pairs[, 1] <- 1:nrow(mp)
    pairs[, 2] <- pi

    if (threshold < min(mp)) {
      stop(paste0("Error: `threshold` is too small for this Matrix Profile. Min: ", round(min(data), 2), ", Max: ", round(max(data), 2)), call. = FALSE)
    }

    # remove excess of arcs
    exclusion.zone <- floor(window.size * exclusion.zone)
    edge.limit <- floor(window.size * edge.limit)
    data[1:edge.limit, ] <- Inf
    data[(data.size - edge.limit + 1):data.size, ] <- Inf

    ind <- which(data < threshold)
    pairs <- pairs[ind, ]

    pairdiff <- pairs[, 1] - pairs[, 2]
    ind <- which(abs(pairdiff) > exclusion.zone)
    pairs <- pairs[ind, ]
  } else {
    data.size <- max(pairs)
  }

  segments <- quality

  z.seq <- seq(0, base::pi, length.out = segments)
  xlim <- c(0, data.size + 1)
  ylim <- c(0, data.size / 2)

  if (is.null(alpha)) {
    alpha <- min(0.5, max(10 / nrow(pairs), 0.03))
  }

  arccolr <- adjustcolor(col, alpha.f = alpha)
  if (length(col) > 1) {
    arccoll <- adjustcolor(col[2], alpha.f = alpha)
  } else {
    arccoll <- adjustcolor(col, alpha.f = alpha)
  }

  # blank plot
  plot(0.5, 0.5,
    type = "n", main = main, xlab = xlab, ylab = ylab,
    xlim = xlim, ylim = ylim, yaxt = "n", ...
  )

  for (i in 1:nrow(pairs)) {
    if (pairs[i, 1] > pairs[i, 2]) {
      arccol <- arccoll
    } else {
      arccol <- arccolr
    }

    x1 <- min(pairs[i, 1], pairs[i, 2])
    x2 <- max(pairs[i, 1], pairs[i, 2])
    center <- (x1 - x2) / 2 + x2
    radius <- (x2 - x1) / 2
    x.seq <- center + radius * cos(z.seq)
    y.seq <- radius * sin(z.seq)
    lines(x.seq, y.seq,
      col = arccol, lwd = lwd, lty = 1, lend = 1
    )
  }

  legend(1, data.size / 2,
    legend = c("Right", "Left"),
    col = adjustcolor(col, alpha.f = 0.5), lty = 1, cex = 0.8, lwd = 5
  )

  invisible(return(pairs))
}
