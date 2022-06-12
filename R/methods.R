#' Model Residuals
#'
#' @description Extract the residual values from an estimated model.
#' @param object an object of class \dQuote{tsvets.estimate}.
#' @param raw raw residuals are the model based values in transformed space 
#' (when either Box Cox or Logistic have been used as transformations).
#' @param ... not currently used.
#' @return An xts matrix of the model residuals.
#' @aliases residuals
#' @method residuals tsvets.estimate
#' @rdname residuals
#' @export
#'
#'
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

#' Model Fitted Values
#'
#' @description Extract the fitted values from an estimated model.
#' @param object an object of class \dQuote{tsvets.estimate}.
#' @param ... not currently used.
#' @aliases fitted
#' @method fitted tsvets.estimate
#' @rdname fitted
#' @export
#'
#'
fitted.tsvets.estimate = function(object, ...)
{
    f <- object$fitted
    colnames(f) <- object$spec$target$y_names
    return(f)
}

#' Extract Model Coefficients
#'
#' @description Extract the estimated coefficients of a model.
#' @param object an object of class \dQuote{tsvets.estimate}.
#' @param ... not currently used.
#' @return a numeric named vector.
#' @aliases coef
#' @method coef tsvets.estimate
#' @rdname coef
#' @export
#'
#'
coef.tsvets.estimate = function(object, ...)
{
    pars <- object$opt$par
    names(pars) <- object$spec$vets_env$parnames
    return(pars)
}

#' Model Decomposition
#'
#' @description Decomposes the estimated model or prediction into its component
#' parts (states).
#' @param object an object of class \dQuote{tsvets.estimate}.
#' @param simplify simplification of the returned states aggregates the level 
#' and slope (if present) into a Trend, Seasonal, X and Irregular components.
#' @param ... not currently used.
#' @return A list of xts matrices with the state components (including error).
#' @aliases tsdecompose
#' @method tsdecompose tsvets.estimate
#' @rdname tsdecompose
#' @export
#'
#'
tsdecompose.tsvets.estimate <- function(object, simplify = FALSE, ...)
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
    if (simplify) {
        Trend <- Level
        if (!is.null(Slope)) Trend <- Trend + Slope
        L <- list(Trend = Trend, Seasonal = Seasonal, X = X, Irregular = residuals(object, raw = TRUE))
    } else {
        L <- list(Level = Level, Slope = Slope, Seasonal = Seasonal, X = X, Irregular = residuals(object, raw = TRUE))
    }
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

#' Model Covariance and Correlation
#'
#' @description Extracts the estimated covariance or correlation matrices.
#' @param object an object of class \dQuote{tsvets.estimate}.
#' @param ... not currently used.
#' @return A Matrix object.
#' @aliases tscov tscor
#' @method tscov tsvets.estimate
#' @rdname correlation
#' @export
#'
#'
tscov.tsvets.estimate = function(object, ...)
{
    # remove any missing values
    E <- object$Error[as.logical(object$spec$target$good_index),,drop = F]
    if (object$spec$dependence$type == "equicorrelation") {
        rho <- object$opt$par[length(object$opt$par)]
        I <- diag(object$spec$vets_env$model[2])
        one <- matrix(1, object$spec$vets_env$model[2], object$spec$vets_env$model[2])
        R <- (1 - rho)*I + rho*one
        v <- diag(apply(E, 2, sd), ncol(E), ncol(E))
        C <-  v %*% R %*% t(v)
    } else if (object$spec$dependence$type == "shrinkage") {
        rho <- object$opt$par[length(object$opt$par)]
        S <- cov(E)
        n <- ncol(E)
        C <- (1 - rho) * S + rho * sum(diag(S))/n * diag(1,n,n)
    } else if (object$spec$dependence$type == "diagonal") {
        # ToDo: need a seperate one for the diagonal
        C <- diag(apply(E, 2, var), ncol(E), ncol(E))
    } else {
        C <- cov(E)
    }
    colnames(C) <- rownames(C) <- object$spec$target$y_names
    return(Matrix(C))
}

#' @method tscor tsvets.estimate
#' @rdname correlation
#' @export
#'
#'
tscor.tsvets.estimate = function(object, ...)
{
    E <- object$Error[as.logical(object$spec$target$good_index),,drop = FALSE]
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

#' Model Log-Likelihood
#'
#' @description Extract the log-likelihood from an estimated model.
#' @param object an object of class \dQuote{tsvets.estimate}.
#' @param ... not currently used.
#' @return Returns an object of class logLik. This is a number with at least one
#' attribute, "df" (degrees of freedom), giving the number of (estimated)
#' parameters in the model.
#' @aliases logLik
#' @method logLik tsvets.estimate
#' @rdname logLik
#' @export
#'
#'
logLik.tsvets.estimate = function(object, ...)
{
    n_states <- ncol(object$States)
    n_pars <- length(object$opt$par)
    # We don't count the covariance for now
    structure(-0.5 * ifelse(is.null(object$opt$negative_llh), NA, object$opt$negative_llh), df = n_states + n_pars, class = "logLik")
}

#' Akaike's An Information Criterion
#'
#' @description Extract the AIC from an estimated model.
#' @param object an object of class \dQuote{tsvets.estimate}.
#' @param ... not currently used.
#' @param k the penalty per parameter to be used; the default k = 2 is the
#' classical AIC.
#' @return a numeric value.
#' @aliases AIC
#' @method AIC tsvets.estimate
#' @rdname AIC
#' @export
#'
#'
AIC.tsvets.estimate <- function(object, ..., k = 2)
{
    LL <- logLik(object)
    np <- attr(LL, "df")
    as.numeric(LL) + k * np
}

#' Performance Metrics
#'
#' @description Performance metrics from an estimated or predicted tsvets model.
#' @param object an object of class \dQuote{tsvets.estimate} or \dQuote{tsvets.predict}
#' @param actual the actual data matched to the dates of the forecasts.
#' @param weights a vector of series importance weights (which sum to 1), used 
#' in generating weighted performance metrics.
#' @param alpha the coverage level for distributional forecast metrics.
#' @param ... not currently used.
#' @aliases tsmetrics
#' @method tsmetrics tsvets.predict
#' @rdname tsmetrics
#' @export
#'
#'
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

#' @method tsmetrics tsvets.estimate
#' @rdname tsmetrics
#' @export
#'
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

#' Model Estimation Summary
#'
#' @description Summary method for class \dQuote{tsvets.estimate}
#' @param object an object of class \dQuote{tsvets.estimate}.
#' @param weights a vector of weights (summing to 1) which is passed to the 
#' \code{tsmetrics} method for generating weighted performance metrics.
#' @param ... not currently used.
#' @return A printout of the parameter summary, model type and some model metrics.
#' @aliases summary
#' @method summary tsvets.estimate
#' @rdname summary
#' @export
#'
#'
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
