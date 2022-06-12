#' Walk Forward Model Backtest
#'
#' @description Generates an expanding window walk forward backtest.
#' @param object an object of class \dQuote{tsvets.spec}.
#' @param start numeric data index from which to start the backtest.
#' @param end numeric data index on which to end the backtest. The backtest will
#' end 1 period before that date in order to have at least 1 out of sample value
#' to compare against.
#' @param h forecast horizon. As the expanding window approaches the \dQuote{end},
#' the horizon will automatically shrink to the number of available out of sample
#' periods.
#' @param alpha optional numeric vector of coverage rates for which to calculate
#' the quantiles.
#' @param aggregate whether to create an aggregate series from the forecasts 
#' (see \code{\link{tsaggregate}} for details).
#' @param weights numeric weights vector for the aggregation.
#' @param solver solver to use.
#' @param autodiff whether to use automatic differentiation for estimation.
#' This makes use of the tsvetsad package.
#' @param trace whether to show the progress bar. The user is expected to have
#' set up appropriate handlers for this using the \dQuote{progressr} package.
#' @param ... additional arguments passed to the \dQuote{auto_clean} function.
#' @return A list with the following data.tables:
#' \itemize{
#' \item prediction : the backtest table with forecasts and actuals by series
#' \item metrics: a summary performance table showing metrics by 
#' forecast horizon and series
#' }
#' @note The function can use parallel functionality as long as the user has
#' set up a \code{\link[future]{plan}} using the future package.
#' @aliases tsbacktest
#' @method tsbacktest tsvets.spec
#' @rdname tsbacktest
#' @export
#'
#'
tsbacktest.tsvets.spec <- function(object, start = floor(NROW(object$target$y_orig)/2), end = NROW(object$target$y_orig),
                                  h = 1, alpha = NULL, aggregate = FALSE, weights = NULL, solver = "nlminb", 
                                  trace = FALSE, autodiff = FALSE, ...)
{
    data <- xts(object$target$y_orig, object$target$index)
    lambda <- sapply(object$transform, function(x) x$lambda)
    frequency <- object$target$frequency
    if (object$xreg$include_xreg) {
        use_xreg <- TRUE
        xreg <- xts(object$xreg$xreg, object$target$index)
    } else {
        use_xreg <- FALSE
        xreg <- NULL
    }
    start_date <- index(data)[start]
    end_date <- index(data)[end - 1]
    seqdates <- index(data[paste0(start_date,"/", end_date)])
    elapsed_time <- function(idx, end_date, start_date) {
        min(h, which(end_date == idx) - which(start_date == idx))
    }
    if (aggregate) {
        if (is.null(aggregate)) {
            stop("\nweights cannot be NULL when aggregate is TRUE")
        }
        weights <- as.numeric(weights)
        if (length(weights) != ncol(data)) stop("\nweights must be a vector of
                                                length equal to ncol data.")
    }
    
    # setup backtest indices
    horizon <- sapply(1:length(seqdates), function(i){
        min(h, elapsed_time(index(data), index(data)[end], seqdates[i]))
    })
    if (!is.null(alpha)) {
        if (any(alpha <= 0)) {
            stop("\nalpha must be strictly positive")
        }
        if (any(alpha >= 1)) {
            stop("\nalpha must be less than 1")
        }
        quantiles <- as.vector(sapply(1:length(alpha), function(k) c(alpha[k]/2, 1 - alpha[k]/2)))
    } else {
        quantiles <- NULL
    }
    if (trace) {
        prog_trace <- progressor(length(seqdates))
    }
    b %<-% future_lapply(1:length(seqdates), function(i) {
        if (trace) prog_trace()
        y_train <- data[paste0("/", seqdates[i])]
        ix <- which(index(data) == seqdates[i])
        y_test <- data[(ix + 1):(ix + horizon[i])]
        if (use_xreg) {
            xreg_train <- xreg[index(y_train)]
            xreg_test <- xreg[index(y_test)]
            xreg_g <- object$xreg$xreg_include
        } else {
            xreg_train <- NULL
            xreg_test <- NULL
            xreg_g <- NULL
        }
        spec <- vets_modelspec(y_train, level = object$model$level, slope = object$model$slope, 
                               damped = object$model$damped, seasonal = object$model$seasonal, 
                               xreg = xreg_train, xreg_include = xreg_g, group = object$model$group,
                               frequency = object$target$frequency, lambda = lambda, lower = 0, 
                               upper = 1, dependence = object$dependence$type)
        mod <- estimate(spec, solver = solver, autodiff = autodiff, ...)
        p <- predict(mod, h = horizon[i], newxreg = xreg_test, forc_dates = index(y_test))
        if (aggregate) {
            ap <- tsaggregate(p, weights = weights)
            dp <- data.table(series = "Aggregate", Level = list(), Slope = list(), Seasonal = list(), X = list(), Error = list(), Predicted = list(ap))
            p$prediction_table <- rbind(p$prediction_table, dp)
            y_test <- cbind(y_test, xts(coredata(y_test) %*% weights, index(y_test)))
            colnames(y_test)[ncol(y_test)] <- "Aggregate"
        }
        if (!is.null(quantiles)) {
            qp <- lapply(1:nrow(p$prediction_table), function(j){
                qpout <- apply(p$prediction_table[j]$Predicted[[1]]$distribution, 2, quantile, quantiles)
                if (length(quantiles) == 1) {
                    qpout <- matrix(qpout, ncol = 1)
                } else{
                    qpout <- t(qpout)
                }
                colnames(qpout) <- paste0("P", round(quantiles*100))
                qpout <- as.data.table(qpout)
                return(qpout)
            })
          qp <- rbindlist(qp)
        }
        out <- lapply(1:nrow(p$prediction_table), function(j){
            data.table(series = p$prediction_table[j]$series, 
                       "estimation_date" = rep(seqdates[i], horizon[i]), 
                       "horizon" = 1:horizon[i], 
                       "size" = rep(nrow(y_train), horizon[i]),
                       "forecast_dates" = as.character(index(y_test)),
                       "forecast" = as.numeric(colMeans(p$prediction_table[j]$Predicted[[1]]$distribution)), 
                       "actual" = as.numeric(y_test[,j]))
        })
        out <- rbindlist(out)
        if (!is.null(quantiles)) out <- cbind(out, qp)
        return(out)
    }, future.packages = c("tsmethods","tsaux","xts","tsvets","data.table"), future.seed = TRUE) 
    
    b <- eval(b)
    b <- rbindlist(b)
    actual <- NULL
    forecast <- NULL
    metrics <- b[,list(MAPE = mape(actual, forecast), 
                       MSLRE = mslre(actual, forecast), 
                       BIAS = bias(actual, forecast),
                       n = .N), by = c("series","horizon")]
    if (!is.null(alpha)) {
        q_names <- matrix(paste0("P", round(quantiles*100)), ncol = 2, byrow = TRUE)
        q <- do.call(cbind, lapply(1:length(alpha), function(i){
            b[,list(mis = mis(actual, get(q_names[i,1]), get(q_names[i,2]), alpha[i])), by = c("series","horizon")]
        }))
        q <- q[,which(grepl("mis",colnames(q))), with = FALSE]
        colnames(q) <- paste0("MIS[",alpha,"]")
        metrics <- cbind(metrics, q)
    }
    return(list(prediction = b, metrics = metrics))
}