residuals.tsvets.estimate = function(object, raw = FALSE, ...)
{
    if (raw) {
        if (!is.null(object$spec$transform[[1]]$lambda)) {
            actual <- object$spec$target$y
            fitted <- do.call(cbind, lapply(1:length(object$spec$transform), function(i)
                object$spec$transform[[i]]$transform(coredata(object$fitted)[,i], object$spec$transform[[i]]$lambda))) 
            r <- actual - fitted
        } else {
            r <- object$spec$target$y_orig - coredata(object$fitted)
        }
    } else {
        r <- object$spec$target$y_orig - coredata(object$fitted)
    }
    r <- xts(r, object$spec$target$index)
    colnames(r) <- object$spec$target$y_names
    return(r)
}

fitted.tsvets.estimate = function(object, ...)
{
    f <- object$fitted
    colnames(f) <- object$spec$target$y_names
    return(f)
}

coef.tsvets.estimate = function(object, ...)
{
    pars <- object$opt$par
    names(pars) <- object$spec$vets_env$parnames
    return(pars)
}

tsdecompose.tsvets.estimate <- function(object, ...)
{
    Level <- xts(object$States[-1, 1:object$spec$vets_env$model[2]], object$spec$target$index)
    colnames(Level) <- paste0("Level[",object$spec$target$y_names,"]")
    k <- object$spec$vets_env$model[2]
    if (object$spec$model$slope != "none") {
        Slope <- xts(object$States[-1, (k + 1):(k + object$spec$vets_env$model[2])], object$spec$target$index)
        colnames(Slope) <- paste0("Slope[",object$spec$target$y_names,"]")
        
        k <- (k + object$spec$vets_env$model[2])
    } else {
        Slope <- NULL
    }
    if (object$spec$model$seasonal != "none") {
        # extend the seasonal index to include the last season which is included in the filter
        Seasonal <- object$States[-1,  -(1:k)]
        Seasonal <- Seasonal[,(ncol(Seasonal) - object$spec$vets_env$model[2] + 1):ncol(Seasonal)]
        Seasonal <- xts(Seasonal, object$spec$target$index)
        colnames(Seasonal) <- paste0("Seasonal[",object$spec$target$y_names,"]")
    } else{
        Seasonal <- NULL
    }
    if (object$spec$xreg$include_xreg) {
        beta <- object$spec$vets_env$Amat[,object$spec$vets_env$X_index[1]:object$spec$vets_env$X_index[2], drop = FALSE]
        xreg <- object$spec$xreg$xreg
        # contribution of weighted X to each Y
        X <- xts(t(beta %*% t(xreg)), object$spec$target$index)
        colnames(X) <- paste0(object$spec$target$y_names,"[X]")
    } else {
        X <- NULL
    }
    L <- list(Level = Level, Slope = Slope, Seasonal = Seasonal, X = X)
    return(L)
}


.tsdecompose_predict <- function(x, i, object, forc_dates, newxreg = NULL, ...)
{
    date_class <- attr(object$spec$target$sampling, "date_class")
    tsd <- tsdecompose(object)
    n <- NCOL(object$spec$target$y_orig)
    Level <- x[,,i]
    if (NCOL(Level) == 1) Level = matrix(Level, ncol = 1)
    colnames(Level) <- forc_dates
    class(Level) <- "tsmodel.distribution"
    attr(Level, "date_class") <- date_class
    k <- n
    Level <- list(original_series = as.zoo(tsd$Level[,i]), distribution = Level)
    class(Level) <- "tsmodel.predict"
    if (object$spec$model$slope != "none") {
        Slope <- x[, ,(k + i)]
        if (NCOL(Slope) == 1) Slope = matrix(Slope, ncol = 1)
        colnames(Slope) <- forc_dates
        class(Slope) <- "tsmodel.distribution"
        attr(Slope, "date_class") <- date_class
        Slope <- list(original_series = as.zoo(tsd$Slope[,i]), distribution = Slope)
        class(Slope) <- "tsmodel.predict"
        k <- k + n
    } else {
        Slope <- NULL
    }
    if (object$spec$model$seasonal != "none") {
        # extend the seasonal index to include the last season which is included in the filter
        s_index <- tail(1:(k + object$spec$target$frequency * n), n)[i]
        Seasonal <- x[, , s_index]
        if (NCOL(Seasonal) == 1) Seasonal = matrix(Seasonal, ncol = 1)
        colnames(Seasonal) <- forc_dates
        class(Seasonal) <- "tsmodel.distribution"
        attr(Seasonal, "date_class") <- date_class
        Seasonal <- list(original_series = as.zoo(tsd$Seasonal[,i]), distribution = Seasonal)
        class(Seasonal) <- "tsmodel.predict"
    } else{
        Seasonal <- NULL
    }
    if (object$spec$xreg$include_xreg) {
        beta <- object$spec$vets_env$Amat[,object$spec$vets_env$X_index[1]:object$spec$vets_env$X_index[2], drop = FALSE]
        X <- beta[i,] %*% newxreg
        if (NCOL(X) == 1) X = matrix(X, ncol = 1)
        colnames(X) <- forc_dates
        class(X) <- "tsmodel.distribution"
        attr(X, "date_class") <- date_class
        X <- list(original_series = as.zoo(tsd$X[,i]), distribution = X)
        class(X) <- "tsmodel.predict"
    } else {
        X <- NULL
    }
    L <- list(Level = Level, Slope = Slope, Seasonal = Seasonal, X = X)
    return(L)
}



.tsdecompose_simulate <- function(x, i, object, sim_dates, ...)
{
    n <- NCOL(object$spec$target$y_orig)
    Level <- x[,,i]
    if (NCOL(Level) == 1) Level = matrix(Level, ncol = length(sim_dates))
    colnames(Level) <- sim_dates
    class(Level) <- "tsmodel.distribution"
    Level <- list(original_series = NULL, distribution = Level)
    class(Level) <- "tsmodel.predict"
    k <- n
    if (object$spec$model$slope != "none") {
        Slope <- x[, ,(k + i)]
        if (NCOL(Slope) == 1) Slope = matrix(Slope, ncol = length(sim_dates))
        colnames(Slope) <- sim_dates
        class(Slope) <- "tsmodel.distribution"
        Slope <- list(original_series = NULL, distribution = Slope)
        class(Slope) <- "tsmodel.predict"
        k <- k + n
    } else {
        Slope <- NULL
    }
    if (object$spec$model$seasonal != "none") {
        # extend the seasonal index to include the last season which is included in the filter
        s_index <- tail(1:(k + object$spec$target$frequency * n), n)[i]
        Seasonal <- x[, , s_index]
        if (NCOL(Seasonal) == 1) Seasonal = matrix(Seasonal, ncol = length(sim_dates))
        colnames(Seasonal) <- sim_dates
        class(Seasonal) <- "tsmodel.distribution"
        Seasonal <- list(original_series = NULL, distribution = Seasonal)
        class(Seasonal) <- "tsmodel.predict"
        
    } else{
        Seasonal <- NULL
    }
    L <- list(Level = Level, Slope = Slope, Seasonal = Seasonal)
    return(L)
}


tscov.tsvets.estimate = function(object, ...)
{
    # remove any missing values
    E <- object$Error[as.logical(object$spec$target$good_index),]
    if (object$spec$dependence$type == "equicorrelation") {
        rho <- object$opt$par[length(object$opt$par)]
        I <- diag(object$spec$vets_env$model[2])
        one <- matrix(1, object$spec$vets_env$model[2], object$spec$vets_env$model[2])
        R <- (1 - rho)*I + rho*one
        v <- diag(apply(E, 2, sd))
        C <-  v %*% R %*% t(v)
    } else if (object$spec$dependence$type == "shrinkage") {
        rho <- object$opt$par[length(object$opt$par)]
        S <- cov(E)
        n <- ncol(E)
        C <- (1 - rho) * S + rho * sum(diag(S))/n * diag(1,n,n)
    } else if (object$spec$dependence$type == "diagonal") {
        # ToDo: need a seperate one for the diagonal
        C <- diag(apply(E, 2, var))
    } else {
        C <- cov(E)
    }
    colnames(C) <- rownames(C) <- object$spec$target$y_names
    return(Matrix(C))
}

tscor.tsvets.estimate = function(object, ...)
{
    E <- object$Error[as.logical(object$spec$target$good_index),]
    if (object$spec$dependence$type == "equicorrelation") {
        rho <- object$opt$par[length(object$opt$par)]
        I <- diag(object$spec$vets_env$model[2])
        one <- matrix(1, object$spec$vets_env$model[2], object$spec$vets_env$model[2])
        R <- (1 - rho)*I + rho*one
    } else if (object$spec$dependence$type == "shrinkage") {
        rho <- object$opt$par[length(object$opt$par)]
        S <- cov(E)
        n <- ncol(E)
        C <- (1 - rho) * S + rho * sum(diag(S))/n * diag(1,n,n)
        R <- cov2cor(C)
    } else if (object$spec$dependence$type == "diagonal") {
        R <- diag(1, ncol(E), ncol(E))
    } else {
        R <- cor(E)
    }
    colnames(R) <- rownames(R) <- object$spec$target$y_names
    return(Matrix(R))
}

logLik.tsvets.estimate = function(object, ...)
{
    n_states <- ncol(object$States)
    n_pars <- length(object$opt$par)
    # We don't count the covariance for now
    structure(-0.5 * ifelse(is.null(object$opt$negative_llh), NA, object$opt$negative_llh), df = n_states + n_pars, class = "logLik")
}

AIC.tsvets.estimate <- function(object, ..., k = 2)
{
    LL <- logLik(object)
    np <- attr(LL, "df")
    as.numeric(LL) + k * np
}

tsmetrics.tsvets.predict = function(object, actual, weights = NULL, alpha = 0.1, ...)
{
    horizon <- ncol(object$prediction_table$Predicted[[1]]$distribution)
    if (nrow(actual) !=  horizon) stop("\nrows of actual not equal to forecast horizon!")
    if (ncol(actual) !=  NROW(object$prediction_table)) stop("\no. of cols of actual not equal to series forecasted!")
    series <- object$prediction_table$series
    if (!is.null(colnames(actual))) { 
        actual <- actual[,series]
    } else {
        colnames(actual) <- series
    }
    if (is.null(weights)) {
        weights <- rep(1/nrow(object$prediction_table), nrow(object$prediction_table))
    } else {
        weights <- as.numeric(weights)
    }
    # MASE, MIS, MAPE and MSLRE
    mu <- do.call(cbind, lapply(1:ncol(actual), function(i){
        matrix(colMeans(object$prediction_table[i]$Predicted[[1]]$distribution), nrow = horizon)
    }))
    orig_series <- do.call(cbind, lapply(1:ncol(actual), function(i){
        matrix(object$prediction_table[i]$Predicted[[1]]$original_series, ncol = 1)
    }))
    frequency <- object$spec$target$frequency
    if (!is.null(alpha)) alpha <- alpha[1]
    itable <- lapply(1:ncol(actual), function(i){
        MAPE <- mape(as.numeric(actual[,i]), as.numeric(mu[,i]))
        BIAS <- bias(as.numeric(actual[,i]), as.numeric(mu[,i]))
        MSLRE <- mslre(as.numeric(actual[,i]), as.numeric(mu[,i]))
        MASE <- mase(actual[,i], mu[,i], orig_series[,i], frequency = frequency)
        if (!is.null(alpha)) {
            MIS <- mis(actual[,i], lower = apply(object$prediction_table[i]$Predicted[[1]]$distribution, 2, quantile, alpha/2), 
                upper = apply(object$prediction_table[i]$Predicted[[1]]$distribution, 2, quantile, 1 - alpha/2), alpha = alpha)
        } else {
            MIS <- NA
        }
        data.table(series = series[i], MASE = MASE, MAPE = MAPE, BIAS = BIAS, MSLRE = MSLRE, MIS = MIS)
    })
    itable <- rbindlist(itable)
    if (ncol(actual) > 1) {
        mtable <- data.table("WAPE" = wape(actual = actual, predicted = mu, weights = weights),
                             "WSLRE" = wslre(actual = actual, predicted = mu, weights = weights))
    } else {
        mtable <- NULL
    }
    return(list(individual = itable, weighted = mtable))
}

tsmetrics.tsvets.estimate = function(object, weights = NULL, ...)
{
    f <- fitted(object)
    y <- object$spec$target$y_orig
    MAPE <- sapply(1:ncol(y), function(i) mape(as.numeric(y[,i]), as.numeric(f[,i])))
    MSLRE <- sapply(1:ncol(y), function(i) mslre(as.numeric(y[,i]), as.numeric(f[,i])))
    MSE <- sapply(1:ncol(y), function(i) mean((as.numeric(y[,i]) - as.numeric(f[,i]))^2))
    mean_MAPE <- mean(MAPE)
    mean_MSLRE <- mean(MSLRE)
    if (is.null(weights)) {
        w <- rep(1/ncol(y), ncol(y))
    } else{
        w <- weights
    }
    if (NCOL(y) > 1) {
        w_MAPE <- wape(y, f, w)
        w_SLRE <- wslre(y, f, w)
    } else {
        w_MAPE <- mean_MAPE
        w_SLRE <- mean_MSLRE
    }
    LL <- logLik(object)
    np <- attr(LL, "df")
    N <- nrow(f) * ncol(f)
    aic <- as.numeric(LL) + 2 * np
    bic <- as.numeric(LL) + log(N) * np
    aicc <- aic + 2 * np * (np + 1)/(N - np - 1)
    data.frame(N = N, no_pars = np, LogLik = as.numeric(LL), AIC = aic, BIC = bic, AICc = aicc, MAPE = mean_MAPE, MSLRE = mean_MSLRE,
               WAPE = w_MAPE, WSLRE = w_SLRE)
}

summary.tsvets.estimate = function(object, weights = NULL, ...)
{
    
    pars <- object$opt$par
    names(pars) <- object$spec$vets_env$parnames
    cat("\n           Vector ETS            \n")
    cat("-----------------------------------\n")
    cat(paste0("Level : ",object$spec$model$level),"\n")
    cat(paste0("Slope : ",object$spec$model$slope),"\n")
    if (object$spec$model$damped != "none") {
        cat(paste0("Dampening : ",object$spec$model$damped),"\n")
    }
    cat(paste0("Seasonal : ",object$spec$model$seasonal),"\n")
    if (object$spec$dependence$type %in% c("equicorrelation","shrinkage")) {
        cat(paste0("Dependence : ",object$spec$dependence$type,"[rho = ",round(pars['rho'],2),"]"),"\n")
    } else {
        cat(paste0("Dependence : ",object$spec$dependence$type),"\n")
    }
    cat(paste0("No. Series : ", ncol(object$spec$target$y_orig)),"\n")
    cat(paste0("No. TimePoints : ", nrow(object$spec$target$y_orig)),"\n")
    
    model_matrices <- p_matrix(object)
    cat("\nParameter Matrices\n")
    cat("\nLevel Matrix\n")
    show(model_matrices$Level_matrix)
    if (!is.null(model_matrices$Slope_matrix)) {
        cat("\nSlope Matrix\n")
        show(model_matrices$Slope_matrix)
    }
    if (!is.null(model_matrices$Phi_matrix)) {
        cat("\nPhi (Dampening) Matrix\n")
        show(model_matrices$Phi_matrix)
    }
    if (!is.null(model_matrices$Seasonal_matrix)) {
        cat("\nSeasonal Matrix\n")
        show(model_matrices$Seasonal_matrix)
    }
    if (object$spec$xreg$include_xreg) {
        cat("\nRegressor Matrix\n")
        show(model_matrices$X_matrix)
    }
    C <- tscor(object)
    if (NCOL(C) == 1) {
        C <- tscov(object)
        cat("\nCovariance Matrix\n")
        show(C)
    } else {
        cat("\nCorrelation Matrix\n")
        show(C)
    }
    m <- tsmetrics(object, weights = weights)
    if (!is.na(m$AIC)) {
        cat("\nInformation Criteria\n")
        print(data.frame(AIC = as.numeric(sprintf(m$AIC, fmt = "%.2f")),
                         BIC = as.numeric(sprintf(m$BIC, fmt = "%.2f")),
                         AICc = as.numeric(sprintf(m$AICc, fmt = "%.2f")), check.names = FALSE),
              row.names = FALSE, right = FALSE)
    }
    cat("\nAccuracy Criteria\n")
    x <- as.data.frame(matrix(c(sprintf(m$MAPE, fmt = "%.4f"), sprintf(m$WAPE, fmt = "%.4f"),
                                sprintf(m$MSLRE, fmt = "%.4f"), sprintf(m$WSLRE, fmt = "%.4f")), 2, 2, byrow = TRUE))
    colnames(x) <- c("Mean","Weighted")
    rownames(x) <- c("MAPE","MSLRE")
    print(x)
}

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


tsaggregate.tsvets.estimate <- function(object, weights = NULL, return_model = FALSE, ...)
{
    lambda <- sapply(object$spec$transform, function(x) x$lambda)
    condition_transform <- is.null(object$spec$transform[[1]]) | (all(lambda == 0) | all(lambda == 1))
    condition_homogeneous <- object$spec$model$level == "common" & object$spec$model$slope %in% c("none","common") & object$spec$model$seasonal %in% c("none","common") & object$spec$model$damped %in% c("none","common")
    condition_model_return <- all(condition_transform, condition_homogeneous)
    # condition guarantees that return_model can be used
    if (return_model & !condition_model_return) {
            stop("\nmodel based aggregation only available for common components (homogeneous model) AND common lambda of NULL, 0 or 1.")
    }
    n <- NCOL(object$spec$target$y)
    if (is.null(weights)) {
        weights <- matrix(1, ncol = n, nrow = 1)
    } else {
        if (length(as.numeric(weights)) != n) stop("\nweights should be a vector of length equal to ncol y.")   
        weights <- matrix(as.numeric(weights), ncol = 1, nrow = n)
    }
    if (!is.null(object$spec$transform[[1]])) {
        if (all(lambda == 0)) {
            actual <- exp(object$spec$target$y %*% weights)
        } else {
            actual <- object$spec$target$y_orig %*% weights
        }
    } else {
        actual <- object$spec$target$y_orig %*% weights
    }
    actual <- xts(actual, object$spec$target$index)
    colnames(actual) <- "aggregate"
    if (return_model) {
        # lambda either NULL (1) or c(0, 1)
        if (object$spec$xreg$include_xreg) {
            xreg <- xts(object$spec$xreg$xreg, object$spec$target$index)
            xreg_include <- matrix(1, ncol = ncol(xreg), nrow = 1)
        } else {
            xreg <- NULL
            xreg_include <- NULL
        }
        spec <- vets_modelspec(actual, level = "common", slope = object$spec$model$slope, damped = object$spec$model$damped, 
                               seasonal = object$spec$model$seasonal, frequency = object$spec$target$frequency, dependence = "diagonal", 
                               lambda = object$spec$transform$lambda[1], xreg = xreg, xreg_include = xreg_include)
        xseed <- object$spec$vets_env$init_states %*% weights
        spec$vets_env$States[,1] <- xseed[1:NROW(spec$vets_env$States)]
        # calculate the variance
        V <- t(weights) %*% as.matrix(tscov(object)) %*% weights
        # reconstruct the estimate object
        pars <- coef(object)
        # weight any beta from regressors
        if (object$spec$xreg$include_xreg) {
            beta <- p_matrix(object)$X_matrix
            beta_index <- which(grepl("^X",names(pars)))
            pars <- pars[-beta_index]
            beta <- as.matrix(t(beta)) %*% weights
            names(beta) <- paste0("X[", colnames(xreg),"][","aggregate","]")
            pars <- c(pars, beta)
        }
        filt <- vets_filter(pars, spec$vets_env)
        if (!is.null(spec$transform$lambda)) {
            fit <- spec$transform$inverse(filt$fitted[-1,,drop = FALSE], spec$transform$lambda)
            fit <- xts(fit, spec$target$index)
        } else {
            fit <- xts(filt$fitted[-1,], spec$target$index)
        }
        spec$vets_env$Amat <- filt$Amat
        spec$vets_env$Fmat <- filt$Fmat
        spec$vets_env$Gmat <- filt$Gmat
        out <- list(fitted = fit, States = filt$States, Error = filt$Error[-1,,drop = FALSE], spec = spec, opt = list(pars = pars, negative_llh = NA), variance = V)
        class(out) <- "tsvets.estimate"
        return(out)
    } else {
        if (!condition_homogeneous) {
            if (!condition_transform) {
                fit <- fitted(object) %*% weights
                out <- cbind(actual, fit)
                colnames(out) <- c("aggregate","fitted")
                return(out)
            } else {
                if (all(object$spec$transform$lambda == 0)) {
                    fit <- exp(log(coredata(object$fitted)) %*% weights)
                } else {
                    fit <- coredata(fitted(object)) %*% weights
                }
                out <- cbind(actual, fit)
                colnames(out) <- c("aggregate","fitted")
                return(out)
            }
        } else {
            if (!condition_transform) {
                fit <- fitted(object) %*% weights
                out <- cbind(actual, fit)
                colnames(out) <- c("aggregate","fitted")
                return(out)
            } else {
                tsd <- tsdecompose(object)
                Level <- coredata(tsd$Level) %*% weights
                a_names <- c("Level")
                if (!is.null(tsd$Slope)) {
                    Slope <- coredata(tsd$Slope) %*% weights
                    a_names <- c(a_names, "Slope")
                } else {
                    Slope <- NULL
                }
                if (!is.null(tsd$Seasonal)) {
                    Seasonal <- coredata(tsd$Seasonal) %*% weights    
                    a_names <- c(a_names, "Seasonal")
                } else {
                    Seasonal <- NULL
                }
                if (!is.null(tsd$X)) {
                    X <- coredata(tsd$X) %*% weights
                    a_names <- c(a_names, "X")
                } else {
                    X <- NULL
                }
                # if lambda = 0 then we can add in logs and exponentiate back to get a multipicative model (e.g. Revenue = Volume x Price)
                if (!is.null(object$spec$transform$lambda)) {
                    if (all(object$spec$transform$lambda == 0)) {
                        fit <- exp(log(coredata(object$fitted)) %*% weights)
                        Err <- (object$spec$target$y_orig/coredata(fitted(object)) - 1) %*% weights
                    } else {
                        fit <- coredata(fitted(object)) %*% weights
                        Err <- coredata(residuals(object)) %*% weights
                    }
                } else {
                    fit <- coredata(fitted(object)) %*% weights
                    Err <- coredata(residuals(object)) %*% weights
                }
                a_names <- c("aggregate","fitted", a_names, "residuals")
                out <- xts(cbind(actual, fit, Level, Slope, Seasonal, X, Err), object$spec$target$index)
                colnames(out) <- a_names
                return(out)   
            }
        }
    }
}


tsaggregate.tsvets.predict <- function(object, weights = NULL, ...)
{
    n <- NCOL(object$spec$target$y)
    if (is.null(weights)) {
        weights <- matrix(1, ncol = n, nrow = 1)
    } else {
        if (length(as.numeric(weights)) != n) stop("\nweights should be a vector of length equal to ncol y.")   
        weights <- matrix(as.numeric(weights), ncol = 1, nrow = n)
    }
    condition_transform <- is.null(object$spec$transform$lambda) | (all(object$spec$transform$lambda == 0) | all(object$spec$transform$lambda == 1))
    condition_homogeneous <- object$spec$model$level == "common" & object$spec$model$slope %in% c("none","common") & object$spec$model$seasonal %in% c("none","common") & object$spec$model$damped %in% c("none","common")
    
    if (!condition_homogeneous) {
        if (condition_transform) {
            f_dates <- colnames(object$prediction_table[1]$Predicted[[1]]$distribution)
            date_class <- attr(object$prediction_table[1]$Predicted[[1]]$distribution, "date_class")
            array_predicted <- array(unlist(lapply(1:nrow(object$prediction_table), function(i) as.matrix(object$prediction_table[i]$Predicted[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$prediction_table)))
            array_predicted <- aperm(array_predicted, perm = c(1,3,2))
            if (all(object$spec$transform$lambda == 0)) {
                distribution <- apply(array_predicted, 3, function(x) exp(log(x) %*% weights))
                original_series <- xts(exp(log(object$spec$target$y_orig) %*% weights), object$spec$target$index)
            } else {
                distribution <- apply(array_predicted, 3, function(x) x %*% weights)
                original_series <- xts(object$spec$target$y_orig %*% weights, object$spec$target$index)
            }
            colnames(distribution) <- f_dates
            class(distribution) <- "tsmodel.distribution"
            attr(distribution, "date_class") <- date_class
            out <- list(original_series = original_series, distribution = distribution)
            class(out) <- "tsmodel.predict"
            return(out)
        } else {
            f_dates <- colnames(object$prediction_table[1]$Predicted[[1]]$distribution)
            date_class <- attr(object$prediction_table[1]$Predicted[[1]]$distribution, "date_class")
            array_predicted <- array(unlist(lapply(1:nrow(object$prediction_table), function(i) as.matrix(object$prediction_table[i]$Predicted[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$prediction_table)))
            array_predicted <- aperm(array_predicted, perm = c(1,3,2))
            distribution <- apply(array_predicted, 3, function(x) x %*% weights)
            colnames(distribution) <- f_dates
            class(distribution) <- "tsmodel.distribution"
            attr(distribution, "date_class") <- date_class
            original_series <- xts(object$spec$target$y_orig %*% weights, object$spec$target$index)
            out <- list(original_series = original_series, distribution = distribution)
            class(out) <- "tsmodel.predict"
            return(out)
        }
    } else {
        # homogenous model so we can aggregate the states
        if (condition_transform) {
            if (!is.null(object$spec$transform$lambda)) {
                # we can aggregate states for lambda = 0, 1, or NULL
                if (all(object$spec$transform$lambda == 0)) {
                    f_dates <- colnames(object$prediction_table[1]$Predicted[[1]]$distribution)
                    date_class <- attr(object$prediction_table[1]$Predicted[[1]]$distribution, "date_class")
                    array_predicted <- array(unlist(lapply(1:nrow(object$prediction_table), function(i) as.matrix(object$prediction_table[i]$Predicted[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$prediction_table)))
                    array_predicted <- aperm(array_predicted, perm = c(1,3,2))
                    distribution <- apply(array_predicted, 3, function(x) exp(log(x) %*% weights))
                    colnames(distribution) <- f_dates
                    class(distribution) <- "tsmodel.distribution"
                    attr(distribution, "date_class") <- date_class
                    original_series <- xts(exp(log(object$spec$target$y_orig) %*% weights), object$spec$target$index)
                } else if (all(object$spec$transform$lambda == 1)) {
                    f_dates <- colnames(object$prediction_table[1]$Predicted[[1]]$distribution)
                    date_class <- attr(object$prediction_table[1]$Predicted[[1]]$distribution, "date_class")
                    array_predicted <- array(unlist(lapply(1:nrow(object$prediction_table), function(i) as.matrix(object$prediction_table[i]$Predicted[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$prediction_table)))
                    array_predicted <- aperm(array_predicted, perm = c(1,3,2))
                    distribution <- apply(array_predicted, 3, function(x) x %*% weights)
                    colnames(distribution) <- f_dates
                    class(distribution) <- "tsmodel.distribution"
                    attr(distribution, "date_class") <- date_class
                    original_series <- xts(object$spec$target$y_orig %*% weights, object$spec$target$index)
                } else {
                    # there is nothing here
                }
                Level <- array(unlist(lapply(1:nrow(object$prediction_table), function(i) as.matrix(object$prediction_table[i]$Level[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$prediction_table)))
                Level <- aperm(Level, perm = c(1,3,2))
                Level <- apply(Level, 3, function(x) x %*% weights)
                colnames(Level) <- f_dates
                class(Level) <- "tsmodel.distribution"
                
                if (!is.null(object$prediction_table[1]$Slope[[1]]$distribution)) {
                    Slope <- array(unlist(lapply(1:nrow(object$prediction_table), function(i) as.matrix(object$prediction_table[i]$Slope[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$prediction_table)))
                    Slope <- aperm(Slope, perm = c(1,3,2))
                    Slope <- apply(Slope, 3, function(x) x %*% weights)
                    colnames(Slope) <- f_dates
                    class(Slope) <- "tsmodel.distribution"
                    attr(Slope, "date_class") <- date_class
                } else {
                    Slope <- NULL
                }
                if (!is.null(object$prediction_table[1]$Seasonal[[1]]$distribution)) {
                    Seasonal <- array(unlist(lapply(1:nrow(object$prediction_table), function(i) as.matrix(object$prediction_table[i]$Seasonal[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$prediction_table)))
                    Seasonal <- aperm(Seasonal, perm = c(1,3,2))
                    Seasonal <- apply(Seasonal, 3, function(x) x %*% weights)
                    colnames(Seasonal) <- f_dates
                    class(Seasonal) <- "tsmodel.distribution"
                    attr(Seasonal, "date_class") <- date_class
                } else {
                    Seasonal <- NULL
                }
                if (!is.null(object$prediction_table[1]$X[[1]]$distribution)) {
                    X <- array(unlist(lapply(1:nrow(object$prediction_table), function(i) as.matrix(object$prediction_table[i]$X[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$prediction_table)))
                    X <- aperm(X, perm = c(1,3,2))
                    X <- apply(X, 3, function(x) x %*% weights)
                    colnames(X) <- f_dates
                    class(X) <- "tsmodel.distribution"
                    attr(X, "date_class") <- date_class
                } else {
                    X <- NULL
                }
                out <- list(original_series = original_series, distribution = distribution, Level = Level, Slope = Slope, Seasonal = Seasonal, X = X)
                class(out) <- "tsmodel.predict"
                return(out) 
            } else {
                f_dates <- colnames(object$prediction_table[1]$Predicted[[1]]$distribution)
                date_class <- attr(object$prediction_table[1]$Predicted[[1]]$distribution, "date_class")
                array_predicted <- array(unlist(lapply(1:nrow(object$prediction_table), function(i) as.matrix(object$prediction_table[i]$Predicted[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$prediction_table)))
                array_predicted <- aperm(array_predicted, perm = c(1,3,2))
                distribution <- apply(array_predicted, 3, function(x) x %*% weights)
                colnames(distribution) <- f_dates
                class(distribution) <- "tsmodel.distribution"
                attr(distribution, "date_class") <- date_class
                original_series <- xts(object$spec$target$y_orig %*% weights, object$spec$target$index)
                Level <- array(unlist(lapply(1:nrow(object$prediction_table), function(i) as.matrix(object$prediction_table[i]$Level[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$prediction_table)))
                Level <- aperm(Level, perm = c(1,3,2))
                Level <- apply(Level, 3, function(x) x %*% weights)
                colnames(Level) <- f_dates
                class(Level) <- "tsmodel.distribution"
                
                if (!is.null(object$prediction_table[1]$Slope[[1]]$distribution)) {
                    Slope <- array(unlist(lapply(1:nrow(object$prediction_table), function(i) as.matrix(object$prediction_table[i]$Slope[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$prediction_table)))
                    Slope <- aperm(Slope, perm = c(1,3,2))
                    Slope <- apply(Slope, 3, function(x) x %*% weights)
                    colnames(Slope) <- f_dates
                    class(Slope) <- "tsmodel.distribution"
                    attr(Slope, "date_class") <- date_class
                } else {
                    Slope <- NULL
                }
                if (!is.null(object$prediction_table[1]$Seasonal[[1]]$distribution)) {
                    Seasonal <- array(unlist(lapply(1:nrow(object$prediction_table), function(i) as.matrix(object$prediction_table[i]$Seasonal[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$prediction_table)))
                    Seasonal <- aperm(Seasonal, perm = c(1,3,2))
                    Seasonal <- apply(Seasonal, 3, function(x) x %*% weights)
                    colnames(Seasonal) <- f_dates
                    class(Seasonal) <- "tsmodel.distribution"
                    attr(Seasonal, "date_class") <- date_class
                } else {
                    Seasonal <- NULL
                }
                if (!is.null(object$prediction_table[1]$X[[1]]$distribution)) {
                    X <- array(unlist(lapply(1:nrow(object$prediction_table), function(i) as.matrix(object$prediction_table[i]$X[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$prediction_table)))
                    X <- aperm(X, perm = c(1,3,2))
                    X <- apply(X, 3, function(x) x %*% weights)
                    colnames(X) <- f_dates
                    class(X) <- "tsmodel.distribution"
                    attr(X, "date_class") <- date_class
                } else {
                    X <- NULL
                }
                out <- list(original_series = original_series, distribution = distribution, Level = Level, Slope = Slope, Seasonal = Seasonal, X = X)
                class(out) <- "tsmodel.predict"
                return(out)
            }
        } else {
            # homogeneous but lambda is not equal across series so we can't aggregate states
            f_dates <- colnames(object$prediction_table[1]$Predicted[[1]]$distribution)
            date_class <- attr(object$prediction_table[1]$Predicted[[1]]$distribution, "date_class")
            array_predicted <- array(unlist(lapply(1:nrow(object$prediction_table), function(i) as.matrix(object$prediction_table[i]$Predicted[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$prediction_table)))
            array_predicted <- aperm(array_predicted, perm = c(1,3,2))
            distribution <- apply(array_predicted, 3, function(x) x %*% weights)
            colnames(distribution) <- f_dates
            class(distribution) <- "tsmodel.distribution"
            attr(distribution, "date_class") <- date_class
            original_series <- xts(object$spec$target$y_orig %*% weights, object$spec$target$index)
            out <- list(original_series = original_series, distribution = distribution)
            class(out) <- "tsmodel.predict"
            return(out)
        }
    }
}


tsaggregate.tsvets.simulate <- function(object, weights = NULL, ...)
{
    n <- NCOL(object$spec$target$y)
    if (is.null(weights)) {
        weights <- matrix(1, ncol = n, nrow = 1)
    } else {
        if (length(as.numeric(weights)) != n) stop("\nweights should be a vector of length equal to ncol y.")   
        weights <- matrix(as.numeric(weights), ncol = 1, nrow = n)
    }
    condition_transform <- is.null(object$spec$transform$lambda) | (all(object$spec$transform$lambda == 0) | all(object$spec$transform$lambda == 1))
    condition_homogeneous <- object$spec$model$level == "common" & object$spec$model$slope %in% c("none","common") & object$spec$model$seasonal %in% c("none","common") & object$spec$model$damped %in% c("none","common")
    
    if (!condition_homogeneous) {
        if (condition_transform) {
            f_dates <- colnames(object$simulation_table[1]$Simulated[[1]]$distribution)
            date_class <- attr(object$simulation_table[1]$Simulated[[1]]$distribution, "date_class")
            array_predicted <- array(unlist(lapply(1:nrow(object$simulation_table), function(i) as.matrix(object$simulation_table[i]$Simulated[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$simulation_table)))
            array_predicted <- aperm(array_predicted, perm = c(1,3,2))
            if (all(object$spec$transform$lambda == 0)) {
                distribution <- apply(array_predicted, 3, function(x) exp(log(x) %*% weights))
            } else {
                distribution <- apply(array_predicted, 3, function(x) x %*% weights)
            }
            colnames(distribution) <- f_dates
            class(distribution) <- "tsmodel.distribution"
            attr(distribution, "date_class") <- date_class
            out <- list(original_series = NULL, distribution = distribution)
            class(out) <- "tsmodel.predict"
            return(out)
        } else {
            f_dates <- colnames(object$simulation_table[1]$Simulated[[1]]$distribution)
            date_class <- attr(object$simulation_table[1]$Simulated[[1]]$distribution, "date_class")
            array_predicted <- array(unlist(lapply(1:nrow(object$simulation_table), function(i) as.matrix(object$simulation_table[i]$Simulated[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$simulation_table)))
            array_predicted <- aperm(array_predicted, perm = c(1,3,2))
            distribution <- apply(array_predicted, 3, function(x) x %*% weights)
            colnames(distribution) <- f_dates
            class(distribution) <- "tsmodel.distribution"
            attr(distribution, "date_class") <- date_class
            out <- list(original_series = NULL, distribution = distribution)
            class(out) <- "tsmodel.predict"
            return(out)
        }
    } else {
        # homogenous model so we can aggregate the states
        if (condition_transform) {
            if (!is.null(object$spec$transform$lambda)) {
                # we can aggregate states for lambda = 0, 1, or NULL
                if (all(object$spec$transform$lambda == 0)) {
                    f_dates <- colnames(object$simulation_table[1]$Simulated[[1]]$distribution)
                    date_class <- attr(object$simulation_table[1]$Simulated[[1]]$distribution, "date_class")
                    array_predicted <- array(unlist(lapply(1:nrow(object$simulation_table), function(i) as.matrix(object$simulation_table[i]$Simulated[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$simulation_table)))
                    array_predicted <- aperm(array_predicted, perm = c(1,3,2))
                    distribution <- apply(array_predicted, 3, function(x) exp(log(x) %*% weights))
                    colnames(distribution) <- f_dates
                    class(distribution) <- "tsmodel.distribution"
                    attr(distribution, "date_class") <- date_class
                    original_series <- xts(exp(log(object$spec$target$y_orig) %*% weights), object$spec$target$index)
                } else if (all(object$spec$transform$lambda == 1)) {
                    f_dates <- colnames(object$simulation_table[1]$Simulated[[1]]$distribution)
                    date_class <- attr(object$simulation_table[1]$Simulated[[1]]$distribution, "date_class")
                    array_predicted <- array(unlist(lapply(1:nrow(object$simulation_table), function(i) as.matrix(object$simulation_table[i]$Simulated[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$simulation_table)))
                    array_predicted <- aperm(array_predicted, perm = c(1,3,2))
                    distribution <- apply(array_predicted, 3, function(x) x %*% weights)
                    colnames(distribution) <- f_dates
                    class(distribution) <- "tsmodel.distribution"
                    attr(distribution, "date_class") <- date_class
                } else {
                    # there is nothing here
                }
                Level <- array(unlist(lapply(1:nrow(object$simulation_table), function(i) as.matrix(object$simulation_table[i]$Level[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$simulation_table)))
                Level <- aperm(Level, perm = c(1,3,2))
                Level <- apply(Level, 3, function(x) x %*% weights)
                colnames(Level) <- f_dates
                class(Level) <- "tsmodel.distribution"
                
                if (!is.null(object$simulation_table[1]$Slope[[1]]$distribution)) {
                    Slope <- array(unlist(lapply(1:nrow(object$simulation_table), function(i) as.matrix(object$simulation_table[i]$Slope[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$simulation_table)))
                    Slope <- aperm(Slope, perm = c(1,3,2))
                    Slope <- apply(Slope, 3, function(x) x %*% weights)
                    colnames(Slope) <- f_dates
                    class(Slope) <- "tsmodel.distribution"
                    attr(Slope, "date_class") <- date_class
                } else {
                    Slope <- NULL
                }
                if (!is.null(object$simulation_table[1]$Seasonal[[1]]$distribution)) {
                    Seasonal <- array(unlist(lapply(1:nrow(object$simulation_table), function(i) as.matrix(object$simulation_table[i]$Seasonal[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$simulation_table)))
                    Seasonal <- aperm(Seasonal, perm = c(1,3,2))
                    Seasonal <- apply(Seasonal, 3, function(x) x %*% weights)
                    colnames(Seasonal) <- f_dates
                    class(Seasonal) <- "tsmodel.distribution"
                    attr(Seasonal, "date_class") <- date_class
                } else {
                    Seasonal <- NULL
                }
                if (!is.null(object$simulation_table[1]$X[[1]]$distribution)) {
                    X <- array(unlist(lapply(1:nrow(object$simulation_table), function(i) as.matrix(object$simulation_table[i]$X[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$simulation_table)))
                    X <- aperm(X, perm = c(1,3,2))
                    X <- apply(X, 3, function(x) x %*% weights)
                    colnames(X) <- f_dates
                    class(X) <- "tsmodel.distribution"
                    attr(X, "date_class") <- date_class
                } else {
                    X <- NULL
                }
                out <- list(original_series = NULL, distribution = distribution, Level = Level, Slope = Slope, Seasonal = Seasonal, X = X)
                class(out) <- "tsmodel.predict"
                return(out) 
            } else {
                f_dates <- colnames(object$simulation_table[1]$Simulated[[1]]$distribution)
                date_class <- attr(object$simulation_table[1]$Simulated[[1]]$distribution, "date_class")
                array_predicted <- array(unlist(lapply(1:nrow(object$simulation_table), function(i) as.matrix(object$simulation_table[i]$Simulated[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$simulation_table)))
                array_predicted <- aperm(array_predicted, perm = c(1,3,2))
                distribution <- apply(array_predicted, 3, function(x) x %*% weights)
                colnames(distribution) <- f_dates
                class(distribution) <- "tsmodel.distribution"
                attr(distribution, "date_class") <- date_class
                Level <- array(unlist(lapply(1:nrow(object$simulation_table), function(i) as.matrix(object$simulation_table[i]$Level[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$simulation_table)))
                Level <- aperm(Level, perm = c(1,3,2))
                Level <- apply(Level, 3, function(x) x %*% weights)
                colnames(Level) <- f_dates
                class(Level) <- "tsmodel.distribution"
                
                if (!is.null(object$simulation_table[1]$Slope[[1]]$distribution)) {
                    Slope <- array(unlist(lapply(1:nrow(object$simulation_table), function(i) as.matrix(object$simulation_table[i]$Slope[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$simulation_table)))
                    Slope <- aperm(Slope, perm = c(1,3,2))
                    Slope <- apply(Slope, 3, function(x) x %*% weights)
                    colnames(Slope) <- f_dates
                    class(Slope) <- "tsmodel.distribution"
                    attr(Slope, "date_class") <- date_class
                } else {
                    Slope <- NULL
                }
                if (!is.null(object$simulation_table[1]$Seasonal[[1]]$distribution)) {
                    Seasonal <- array(unlist(lapply(1:nrow(object$simulation_table), function(i) as.matrix(object$simulation_table[i]$Seasonal[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$simulation_table)))
                    Seasonal <- aperm(Seasonal, perm = c(1,3,2))
                    Seasonal <- apply(Seasonal, 3, function(x) x %*% weights)
                    colnames(Seasonal) <- f_dates
                    class(Seasonal) <- "tsmodel.distribution"
                    attr(Seasonal, "date_class") <- date_class
                } else {
                    Seasonal <- NULL
                }
                if (!is.null(object$simulation_table[1]$X[[1]]$distribution)) {
                    X <- array(unlist(lapply(1:nrow(object$simulation_table), function(i) as.matrix(object$simulation_table[i]$X[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$simulation_table)))
                    X <- aperm(X, perm = c(1,3,2))
                    X <- apply(X, 3, function(x) x %*% weights)
                    colnames(X) <- f_dates
                    class(X) <- "tsmodel.distribution"
                    attr(X, "date_class") <- date_class
                } else {
                    X <- NULL
                }
                out <- list(original_series = NULL, distribution = distribution, Level = Level, Slope = Slope, Seasonal = Seasonal, X = X)
                class(out) <- "tsmodel.predict"
                return(out)
            }
        } else {
            # homogeneous but lambda is not equal across series so we can't aggregate states
            f_dates <- colnames(object$simulation_table[1]$Simulated[[1]]$distribution)
            date_class <- attr(object$simulation_table[1]$Simulated[[1]]$distribution, "date_class")
            array_predicted <- array(unlist(lapply(1:nrow(object$simulation_table), function(i) as.matrix(object$simulation_table[i]$Simulated[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$simulation_table)))
            array_predicted <- aperm(array_predicted, perm = c(1,3,2))
            distribution <- apply(array_predicted, 3, function(x) x %*% weights)
            colnames(distribution) <- f_dates
            class(distribution) <- "tsmodel.distribution"
            attr(distribution, "date_class") <- date_class
            out <- list(original_series = NULL, distribution = distribution)
            class(out) <- "tsmodel.predict"
            return(out)
        }
    }
}