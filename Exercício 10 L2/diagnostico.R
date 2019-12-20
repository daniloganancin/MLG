diagnostico <- function(object, which = c(1:3, 5), data = NULL, colour = "#444444", 
                        size = NULL, linetype = NULL, alpha = NULL, fill = NULL, 
                        shape = NULL, label = TRUE, label.label = ".label", label.colour = "#000000", 
                        label.alpha = NULL, label.size = NULL, label.angle = NULL, 
                        label.family = NULL, label.fontface = NULL, label.lineheight = NULL, 
                        label.hjust = NULL, label.vjust = NULL, label.repel = FALSE, 
                        label.n = 3, smooth.colour = "#0000FF", smooth.linetype = "solid", 
                        ad.colour = "#888888", ad.linetype = "dashed", ad.size = 0.2, 
                        nrow = NULL, ncol = NULL, xlab=NULL, ylab=NULL, title=NULL, ...) 
{
  p1 <- p2 <- p3 <- p4 <- p5 <- p6 <- NULL
  dropInf <- function(x, h) {
    if (any(isInf <- h >= 1)) {
      warning(gettextf("not plotting observations with leverage one:\n  %s", 
                       paste(which(isInf), collapse = ", ")), call. = FALSE, 
              domain = NA)
      x[isInf] <- NaN
    }
    x
  }
  show <- rep(FALSE, 6)
  show[which] <- TRUE
  if (is.null(data)) {
    plot.data <- ggplot2::fortify(object)
  }
  else {
    plot.data <- ggplot2::fortify(object, data = data)
  }
  n <- nrow(plot.data)
  plot.data$.index <- 1:n
  plot.data$.label <- rownames(plot.data)
  is_glm <- inherits(object, "glm")
  r <- residuals(object)
  w <- weights(object)
  if (any(show[2L:6L])) {
    s <- if (inherits(object, "rlm")) {
      object$s
    }
    else if (is_glm) {
      sqrt(summary(object)$dispersion)
    }
    else {
      sqrt(stats::deviance(object)/stats::df.residual(object))
    }
    hii <- stats::lm.influence(object, do.coef = FALSE)$hat
    if (any(show[2L:3L])) {
      plot.data$.wresid <- if (is.null(w)) {
        r
      }
      else {
        sqrt(w) * r
      }
      plot.data$.wstdresid <- plot.data$.wresid/(s * sqrt(1 - 
                                                            hii))
    }
    if (show[2L]) {
      ylim <- range(plot.data$.wstdresid, na.rm = TRUE)
      ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
      qn <- stats::qqnorm(plot.data$.wstdresid, ylim = ylim, 
                          plot.it = FALSE)
      plot.data$.qqx <- qn$x
      plot.data$.qqy <- qn$y
    }
  }
  label.fitted <- ifelse(is_glm, "Predicted values", "Fitted values")
  label.y23 <- ifelse(is_glm, "Std. deviance resid.", "Standardized residuals")
  if (is.logical(shape) && !shape) {
    if (missing(label)) {
      label <- TRUE
    }
    if (missing(label.n)) {
      label.n <- nrow(plot.data)
    }
  }
  plot.data <- ggfortify:::flatten(plot.data)
  if (label.n > 0L) {
    if (show[1L]) {
      r.data <- dplyr::arrange_(plot.data, "dplyr::desc(abs(.resid))")
      r.data <- utils::head(r.data, label.n)
    }
    if (".wresid" %in% colnames(plot.data)) {
      wr.data <- dplyr::arrange_(plot.data, "dplyr::desc(abs(.wresid))")
      wr.data <- utils::head(wr.data, label.n)
    }
    if (any(show[4L:6L])) {
      cd.data <- dplyr::arrange_(plot.data, "dplyr::desc(abs(.cooksd))")
      cd.data <- utils::head(cd.data, label.n)
    }
  }
  .smooth <- function(x, y) {
    stats::lowess(x, y, f = 2/3, iter = 3)
  }
  .decorate.label <- function(p, data) {
    if (label & label.n > 0) {
      p <- ggfortify:::plot_label(p = p, data = data, label = label, 
                                  label.label = label.label, label.colour = label.colour, 
                                  label.alpha = label.alpha, label.size = label.size, 
                                  label.angle = label.angle, label.family = label.family, 
                                  label.fontface = label.fontface, label.lineheight = label.lineheight, 
                                  label.hjust = label.hjust, label.vjust = label.vjust, 
                                  label.repel = label.repel)
    }
    p
  }
  .decorate.plot <- function(p, xlab = NULL, ylab = NULL, title = NULL) {
    p + ggplot2::xlab(xlab) + ggplot2::ylab(ylab) + ggplot2::ggtitle(title)
  }
  smoother_m <- ggplot2::aes_string(x = "x", y = "y")
  if (show[1L]) {
    t1 <- "Residuals vs Fitted"
    mapping <- ggplot2::aes_string(x = ".fitted", y = ".resid")
    smoother <- .smooth(plot.data$.fitted, plot.data$.resid)
    smoother <- as.data.frame(smoother)
    p1 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
    if (!is.logical(shape) || shape) {
      p1 <- p1 + ggfortify:::geom_factory(geom_point, plot.data, colour = colour, 
                                          size = size, linetype = linetype, alpha = alpha, 
                                          fill = fill, shape = shape) + theme_minimal()
    }
    p1 <- p1 + ggplot2::geom_line(data = smoother, mapping = smoother_m, 
                                  colour = smooth.colour, linetype = smooth.linetype) + 
      ggplot2::geom_hline(yintercept = 0L, linetype = ad.linetype, 
                          size = ad.size, colour = ad.colour)
    p1 <- .decorate.label(p1, r.data)
    if (is.null(xlab)) { xlabel.resfit <- label.fitted
    } else {xlabel.resfit <- xlab$resfit}
    if (is.null(ylab)) { ylabel.resfit <- "Residuals"
    } else {ylabel.resfit <- ylab$resfit}
    if (is.null(title)) { title.resfit <- t1
    } else {title.resfit <- title$resfit}
    p1 <- .decorate.plot(p1, xlab = xlabel.resfit, ylab = ylabel.resfit, 
                         title = title.resfit)
  }
  if (show[2L]) {
    t2 <- "Normal Q-Q"
    qprobs <- c(0.25, 0.75)
    qy <- stats::quantile(plot.data$.wstdresid, probs = qprobs, 
                          names = FALSE, type = 7, na.rm = TRUE)
    qx <- stats::qnorm(qprobs)
    slope <- diff(qy)/diff(qx)
    int <- qy[1L] - slope * qx[1L]
    mapping <- ggplot2::aes_string(x = ".qqx", y = ".qqy")
    p2 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
    if (!is.logical(shape) || shape) {
      p2 <- p2 + ggfortify:::geom_factory(geom_point, plot.data, colour = colour, 
                                          size = size, linetype = linetype, alpha = alpha, 
                                          fill = fill, shape = shape) + theme_minimal()
    }
    p2 <- p2 + ggplot2::geom_abline(intercept = int, slope = slope, 
                                    linetype = ad.linetype, size = ad.size, colour = ad.colour)
    p2 <- .decorate.label(p2, wr.data)
    if (is.null(xlab)) { xlabel.qqplot <- "Theoretical Quantiles"
    } else {xlabel.qqplot <- xlab$qqplot}
    if (is.null(ylab)) { ylabel.qqplot <- label.y23
    } else {ylabel.qqplot <- ylab$qqplot}
    if (is.null(title)) { title.qqplot <- t2
    } else {title.qqplot <- title$qqplot}
    p2 <- .decorate.plot(p2, xlab = xlabel.qqplot, ylab = ylabel.qqplot, 
                         title = title.qqplot)
  }
  if (show[3L]) {
    t3 <- "Scale-Location"
    mapping <- ggplot2::aes_string(x = ".fitted", y = "sqrt(abs(.wstdresid))")
    smoother <- .smooth(plot.data$.fitted, sqrt(abs(plot.data$.wstdresid)))
    smoother <- as.data.frame(smoother)
    p3 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
    if (!is.logical(shape) || shape) {
      p3 <- p3 + ggfortify:::geom_factory(geom_point, plot.data, colour = colour, 
                                          size = size, linetype = linetype, alpha = alpha, 
                                          fill = fill, shape = shape) + theme_minimal()
    }
    p3 <- p3 + ggplot2::geom_line(data = smoother, mapping = smoother_m, 
                                  colour = smooth.colour, linetype = smooth.linetype)
    p3 <- .decorate.label(p3, wr.data)
    label.y3 <- ifelse(is_glm, expression(sqrt(abs(`Std. deviance resid.`))), 
                       expression(sqrt(abs(`Standardized residuals`))))
    if (is.null(xlab)) { xlabel.scaleloc <- label.fitted
    } else {xlabel.scaleloc <- xlab$scaleloc}
    if (is.null(ylab)) { ylabel.scaleloc <- label.y3
    } else {ylabel.scaleloc <- ylab$scaleloc}
    if (is.null(title)) { title.scaleloc <- t3
    } else {title.scaleloc <- title$scaleloc}
    p3 <- .decorate.plot(p3, xlab = xlabel.scaleloc, ylab = ylabel.scaleloc, 
                         title = title.scaleloc) 
  }
  if (show[4L]) {
    t4 <- "Cook's distance"
    mapping <- ggplot2::aes_string(x = ".index", y = ".cooksd", 
                                   ymin = 0, ymax = ".cooksd")
    p4 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
    if (!is.logical(shape) || shape) {
      p4 <- p4 + ggfortify:::geom_factory(geom_linerange, plot.data, 
                                          colour = colour, size = size, linetype = linetype, 
                                          alpha = alpha, fill = fill, shape = shape) + theme_minimal()
    }
    p4 <- .decorate.label(p4, cd.data)
    if (is.null(xlab)) { xlabel.cook <- "Obs. Number"
    } else {xlabel.cook <- xlab$cook}
    if (is.null(ylab)) { ylabel.cook <- "Cook's distance"
    } else {ylabel.cook <- ylab$cook}
    if (is.null(title)) { title.cook <- t4
    } else {title.cook <- title$cook}
    p4 <- .decorate.plot(p4, xlab = xlabel.cook , ylab = ylabel.cook, 
                         title = title.cook)
  }
  if (show[5L]) {
    t5 <- "Residuals vs Leverage"
    mapping <- ggplot2::aes_string(x = ".hat", y = ".stdresid")
    smoother <- .smooth(plot.data$.hat, plot.data$.stdresid)
    smoother <- as.data.frame(smoother)
    p5 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
    if (!is.logical(shape) || shape) {
      p5 <- p5 + ggfortify:::geom_factory(geom_point, plot.data, colour = colour, 
                                          size = size, linetype = linetype, alpha = alpha, 
                                          fill = fill, shape = shape) + theme_minimal()
    }
    p5 <- p5 + ggplot2::geom_line(data = smoother, mapping = smoother_m, 
                                  colour = smooth.colour, linetype = smooth.linetype) + 
      ggplot2::geom_hline(yintercept = 0L, linetype = ad.linetype, 
                          size = ad.size, colour = ad.colour) + ggplot2::expand_limits(x = 0)
    p5 <- .decorate.label(p5, cd.data)
    label.y5 <- ifelse(is_glm, "Std. Pearson resid.", "Standardized Residuals")
    if (is.null(xlab)) { xlabel.reslev <- "Leverage"
    } else {xlabel.reslev <- xlab$reslev}
    if (is.null(ylab)) { ylabel.reslev <- label.y5
    } else {ylabel.reslev <- ylab$reslev}
    if (is.null(title)) { title.reslev <- t5
    } else {title.reslev <- title$reslev}
    p5 <- .decorate.plot(p5, xlab = xlabel.reslev, ylab = ylabel.reslev, 
                         title = title.reslev)
  }
  if (show[6L]) {
    t6 <- "Cook's dist vs Leverage"
    mapping <- ggplot2::aes_string(x = ".hat", y = ".cooksd")
    smoother <- .smooth(plot.data$.hat, plot.data$.cooksd)
    smoother <- as.data.frame(smoother)
    p6 <- ggplot2::ggplot(data = plot.data, mapping = mapping)
    if (!is.logical(shape) || shape) {
      p6 <- p6 + ggfortify:::geom_factory(geom_point, plot.data, colour = colour, 
                                          size = size, linetype = linetype, alpha = alpha, 
                                          fill = fill, shape = shape) + theme_minimal()
    }
    p6 <- p6 + ggplot2::geom_line(data = smoother, mapping = smoother_m, 
                                  colour = smooth.colour, linetype = smooth.linetype) + 
      ggplot2::expand_limits(x = 0, y = 0)
    p6 <- .decorate.label(p6, cd.data)
    if (is.null(xlab)) { xlabel.cooklev <- "Leverage"
    } else {xlabel.cooklev <- xlab$cooklev}
    if (is.null(ylab)) { ylabel.cooklev <- "Cook's distance"
    } else {ylabel.cooklev <- ylab$cooklev}
    if (is.null(title)) { title.cooklev <- t6
    } else {title.cooklev <- title$cooklev}
    p6 <- .decorate.plot(p6, xlab = xlabel.cooklev, ylab = ylabel.cooklev, 
                         title = title.cooklev)
    g <- dropInf(hii/(1 - hii), hii)
    p <- length(stats::coef(object))
    bval <- pretty(sqrt(p * plot.data$.cooksd/g), 5)
    for (i in seq_along(bval)) {
      bi2 <- bval[i]^2
      p6 <- p6 + ggplot2::geom_abline(intercept = 0, slope = bi2, 
                                      linetype = ad.linetype, size = ad.size, colour = ad.colour)
    }
  }
  if (is.null(ncol)) {
    ncol <- 0
  }
  if (is.null(nrow)) {
    nrow <- 0
  }
  plot.list <- list(p1, p2, p3, p4, p5, p6)[which]
  new("ggmultiplot", plots = plot.list, nrow = nrow, ncol = ncol)
}




