tsbacktest.tsvets.spec <- function(object, start = floor(NROW(object$target$y_orig)/2), end = NROW(object$target$y_orig),
                                  h = 1, alpha = NULL, cores = 1, save_output = FALSE, 
                                  save_dir = "~/tmp/", solver = "nlminb", trace = FALSE, ...)
{
    if (save_output) {
        if (is.null(save_dir)) {
            stop("save_dir cannot be NULL when save.output is TRUE")
        }
        if (!dir.exists(save_dir)) {
            stop("save_dir does not exit. Create first and then resubmit")
        }
    }
    data <- xts(object$target$y_orig, object$target$index)
    lambda <- object$transform$lambda
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
    i <- 1
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    clusterExport(cl, "object", envir = environment())
    clusterExport(cl, "data", envir = environment())
    clusterExport(cl, "seqdates", envir = environment())
    clusterExport(cl, "horizon", envir = environment())
    clusterExport(cl, "use_xreg", envir = environment())
    clusterExport(cl, "xreg", envir = environment())
    if (trace) {
        iterations <- length(seqdates)
        pb <- txtProgressBar(max = iterations, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
    } else {
        opts <- NULL
    }
    b <- foreach(i = 1:length(seqdates), .packages = c("tsmethods","tsaux","xts","tsvets","data.table"), .options.snow = opts) %dopar% {
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
                               frequency = object$target$frequency, lambda = NA, lambda_lower = 0, 
                               lambda_upper = 1, dependence = object$dependence$type, cores = 1)
        mod <- estimate(spec, solver = solver)
        p <- predict(mod, h = horizon[i], newxreg = xreg_test, forc_dates = index(y_test))
        if (save_output) {
            saveRDS(mod, file = paste0(save_dir,"/model_", seqdates[i], ".rds"))
            saveRDS(p, file = paste0(save_dir,"/predict_", seqdates[i], ".rds"))
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
    }
    stopCluster(cl)
    if (trace) {
        close(pb)
    }
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