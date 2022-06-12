#' Object Plots
#'
#' @description Plots for objects generated from the tsvets functions.
#' @param x an object of class \dQuote{tsvets.estimate} or \dQuote{tsvets.predict}.
#' @param y not used.
#' @param type the type of plot with a choice of \dQuote{fitted}, \dQuote{states} 
#' and \dQuote{residuals}.
#' @param series the selection of series to plot, with a maximum of 10 at a time 
#' for object of class \dQuote{tsvets.estimate}, and the single series for the 
#' object of class \dQuote{tsvets.predict}
#' @param n_original number of actual data points to include in the prediction 
#' plot.
#' @param ... additional arguments passed to the underlying plot function.
#' @aliases plot
#' @method plot tsvets.estimate
#' @rdname plot
#' @export
#'
#'
plot.tsvets.estimate = function(x, y = NULL, type = c("fitted", "states", "residuals"), series = 1:min(10, ncol(fitted(x))), ...)
{
    # add option for custom colors
    opar <- par()
    opar$cin <- NULL
    opar$cra <- NULL
    opar$csi <- NULL
    opar$cxy <- NULL
    opar$din <- NULL
    opar$page <- NULL
    type <- match.arg(type[1L], choices = c("fitted", "states", "residuals"), several.ok = FALSE)
    n_series <- 1L:ncol(x$fitted)
    series_index <- x$spec$target$index
    if (!all(series %in% n_series)) {
        stop("\nvalues in series do not match actual series")
    }
    if (length(series) > 10) warning("\nMore than 10 series selected. Truncating to the first 10.")
    series <- series[1:min(10,length(series))]
    if (type == "fitted") {
        k <- length(series)
        par(mfrow = n2mfrow(k), mar = c(2.5,2.5,2.5,2.5))
        for (i in 1L:k) {
            plot(zoo(x$spec$target$y_orig[,series[i]], series_index), main = x$spec$target$y_names[i], ylab = "", xlab = "")
            lines(as.zoo(x$fitted[,series[i]]), col = "tomato1", lwd = 0.8)
            grid()
        }
    } else if (type == "states") {
        tsd <- tsdecompose(x)
        k <- 1
        n <- ncol(x$fitted)
        if (x$spec$model$slope != "none") {
            k <- k + 1
        }
        if (x$spec$model$seasonal != "none") {
            k <- k + 1
        }
        if (x$spec$xreg$include_xreg) {
            k <- k + 1
        }
        if (k == 1) {
            m <- matrix(c(1,2), nrow = 2, ncol = 1, byrow = TRUE)
            layout(mat = m, heights = c(0.4,0.2))
        } else if (k == 2) {
            m <- matrix(c(1,1,2,2,3,3), nrow = 3, ncol = 2,byrow = TRUE)
            layout(mat = m, heights = c(0.4,0.4,0.2))
        } else if (k == 3) {
            m <- matrix(c(1,1,2,2,3,3,4,4), nrow = 4, ncol = 2, byrow = TRUE)
            layout(mat = m)
        } else {
            m <- matrix(c(1,1,2,2,3,3,4,4,5,5), nrow = 5, ncol = 2, byrow = TRUE)
            layout(mat = m)
        }
        par(bg = "white", mar = c(0.5,2,0.5,3))
        plot_colors <- viridis_pal(option = "D", end = 0.8)(length(series))
        ylim <- c(min(coredata(tsd$Level[,series])),max(coredata(tsd$Level[,series])))
        if (k == 1) xaxt <- "s" else xaxt <- "n"
        plot(index(tsd$Level), coredata(tsd$Level[,series[1]]), col = plot_colors, type = "l", lty = 1, ylim = ylim, main = "", ylab = "", xlab = "", xaxt = xaxt, cex.axis = 0.8)
        for (i in 2:length(series)) lines(index(tsd$Level), coredata(tsd$Level[,series[i]]), col = plot_colors[i])
        mtext("Level", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
        grid()
        if (x$spec$model$slope != "none") {
            if (k == 2) xaxt <- "s" else xaxt <- "n"
            par(bg = "white", mar = c(0.5,2,0.1,3))
            ylim <- c(min(coredata(tsd$Slope[,series])),max(coredata(tsd$Slope[,series])))
            plot(index(tsd$Slope), coredata(tsd$Slope[,series[1]]), col = plot_colors, type = "l", lty = 1, ylim = ylim, main = "", ylab = "", xlab = "", xaxt = xaxt, cex.axis = 0.8)
            for (i in 2:length(series)) lines(index(tsd$Slope), coredata(tsd$Slope[,series[i]]), col = plot_colors[i])
            mtext("Slope", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
            grid()
        }
        if (x$spec$model$seasonal != "none") {
            if (x$spec$xreg$include_xreg) xaxt <- "n" else xaxt <- "s"
            par(bg = "white", mar = c(0.5,2,0.1,3))
            ylim <- c(min(coredata(tsd$Seasonal[,series])),max(coredata(tsd$Seasonal[,series])))
            plot(index(tsd$Seasonal), coredata(tsd$Seasonal[,series[1]]), col = plot_colors, type = "l", lty = 1, ylim = ylim, main = "", ylab = "", xlab = "", xaxt = xaxt, cex.axis = 0.8)
            for (i in 2:length(series)) lines(index(tsd$Seasonal), coredata(tsd$Seasonal[,series[i]]), col = plot_colors[i])
            mtext("Seasonal", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
            grid()
        }
        if (x$spec$xreg$include_xreg) {
            xaxt <- "s"
            par(bg = "white", mar = c(0.5,2,0.1,3))
            ylim <- c(min(coredata(tsd$X[,series])),max(coredata(tsd$X[,series])))
            plot(index(tsd$X), coredata(tsd$X[,series[1]]), col = plot_colors, type = "l", lty = 1, ylim = ylim, main = "", ylab = "", xlab = "", xaxt = xaxt, cex.axis = 0.8)
            for (i in 2:length(series)) lines(index(tsd$X), coredata(tsd$X[,series[i]]), col = plot_colors[i])
            mtext("Regressors", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
            grid()
        }
        par(bg = "white", mar = c(2,2,3,3))
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        legend(x = "top", inset = 0, legend = x$spec$target$y_names[series], col = plot_colors, lwd = 5, cex = 0.7, ncol = ifelse(n > 5, 2, 1))
    } else if (type == "residuals") {
        k <- length(series)
        par(mfrow = n2mfrow(k), mar = c(2.5,2.5,2.5,2.5))
        for (i in 1L:k) {
            plot(zoo(x$Error[,series[i]], series_index), main = x$spec$target$y_names[series[i]], ylab = "", xlab = "")
            grid()
        }
    }
    suppressWarnings(par(opar))
}

#' @method plot tsvets.predict
#' @rdname plot
#' @export
plot.tsvets.predict = function(x, y = NULL, series = 1, n_original = NULL, ...)
{
    # add option for custom colors
    opar <- par()
    opar$cin <- NULL
    opar$cra <- NULL
    opar$csi <- NULL
    opar$cxy <- NULL
    opar$din <- NULL
    opar$page <- NULL
    k <- 2
    state_names <- c("Level")
    if (!is.null(x$prediction_table$Slope[[1]][[1]])) {
        k <- k + 1
        state_names <- c(state_names, "Slope")
    }
    if (!is.null(x$prediction_table$Seasonal[[1]][[1]])) {
        k <- k + 1
        state_names <- c(state_names, "Seasonal")
    }
    if (!is.null(x$prediction_table$X[[1]][[1]])) {
        k <- k + 1
        state_names <- c(state_names, "X")
    }
    n_series <- 1L:nrow(x$prediction_table)
    if (is.null(series)) series <- 1
    if (length(series) > 1) {
        series <- series[1]
        warning("\nseries must be of length 1. Truncating to the first value.")
    }
    if (!all(series %in% n_series)) {
        stop("\nvalues in series do not match actual series")
    }
    select <- x$prediction_table[series]
    n <- k
    m <- matrix(sort(rep(1:n, times = 2)), nrow = n, ncol = 2, byrow = TRUE)
    plot_colors <- c("steelblue","steelblue","steelblue","steelblue","grey")
    layout(mat = m)
    par(bg = "white", mar = c(0.5,2,2,3))
    plot(select$Predicted[[1]], n_original = n_original, gradient_color = plot_colors[1], interval_color = plot_colors[1], x_axes = FALSE, main = select$series)
    mtext("Predicted", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
    par(bg = "white", mar = c(0.5,2,0.5,3))
    for (i in 1:length(state_names)) {
        if (i == length(state_names)) {
            xax <- TRUE
            par(bg = "white", mar = c(2, 2,0.5,3))
        } else {
            xax <- FALSE
        }
        plot(select[,state_names[i], with = FALSE][[1]][[1]], n_original = n_original, gradient_color = plot_colors[i + 1], interval_color = plot_colors[i + 1], x_axes = xax)
        mtext(state_names[i], side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
    }
    suppressWarnings(par(opar))
}