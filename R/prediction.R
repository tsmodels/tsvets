predict.tsvets.estimate <- function(object, h = 12, newxreg = NULL, nsim = 1000, forc_dates = NULL, init_states = NULL, ...)
{
    if (!is.null(forc_dates)) {
        if (h != length(forc_dates))
            stop("\nforc_dates must have length equal to h")
    }
    if (!object$spec$xreg$include_xreg) {
        X <- matrix(0, ncol = 1, nrow = h + 1)
        if (is.null(forc_dates)) {
            forc_dates = future_dates(tail(object$spec$target$index, 1), frequency = object$spec$target$sampling, n = h)
        }
    } else {
        if (!is.null(newxreg)) {
            if (is.xts(newxreg)) {
                forc_dates = index(newxreg)
            } else {
                forc_dates = future_dates(tail(object$spec$target$index, 1), frequency = object$spec$target$sampling, n = h)
            }
            h <- NROW(newxreg)
            X <- rbind(matrix(0, ncol = ncol(newxreg), nrow = 1), coredata(newxreg))
        } else {
            if (is.null(forc_dates)) {
                forc_dates = future_dates(tail(object$spec$target$index, 1), frequency = object$spec$target$sampling, n = h)
            }
            warning("\nxreg use in estimation but newxreg is NULL...setting to zero")
            X <- matrix(0, ncol = ncol(object$spec$xreg$xreg), nrow = h + 1)
        }
    }
    X <- t(X)
    S <- as.matrix(tscov(object))
    if (!is.positive.definite(S)) {
        warning("\ncovariance is not positive definite...making positive definite.")
        S <- make.positive.definite(S)
    }
    date_class <- attr(object$spec$target$sampling, "date_class")
    xseed <- tail(object$States, 1)
    if (!is.null(init_states)) {
        if (length(as.vector(init_states) != ncol(xseed))) {
            stop(paste0("\ninit_states must be a vector of length ", ncol(xseed)))
        } else {
            xseed <- matrix(as.numeric(init_states), nrow = 1, ncol = ncol(xseed))
        }
    }
    xseed <- t(xseed)
    Amat <- object$spec$vets_env$Amat
    Amat[object$spec$vets_env$matrix_index] <- object$opt$par[object$spec$vets_env$parameter_index]
    Fmat <- object$spec$vets_env$Fmat
    if (object$spec$vets_env$Phi_index[1] > 0) {
        Fmat[which(as.logical(is.na(Fmat)))] <- Amat[,object$spec$vets_env$Phi_index[1]:object$spec$vets_env$Phi_index[2]]
    }
    if (object$spec$vets_env$X_index[1] >= 0) {
        beta <- Amat[,object$spec$vets_env$X_index[1]:object$spec$vets_env$X_index[2], drop = FALSE]
    } else {
        beta <- matrix(1, ncol = object$spec$vets_env$model[2], nrow = 1)
    }
    Amat <- t(Amat[,object$spec$vets_env$Amat_index[1]:object$spec$vets_env$Amat_index[2], drop = FALSE])
    
    Gmat <- object$spec$vets_env$Gmat
    Hmat <- object$spec$vets_env$Hmat
    model <- object$spec$vets_env$model
    ##################################################################
    model <- c(h, nsim, NCOL(object$fitted), ifelse(object$spec$xreg$include_xreg, 1, 0))
    f <- vets_cpp_predict(model, S, Amat, Fmat, Hmat, Gmat, istate = xseed, X = X, beta)
    Y <- aperm(f$Y, map = list(3,2,1))
    E <- aperm(f$Error, map = list(3,2,1))
    States <- aperm(f$States, list(3,2,1))
    forc_dates <- as.character(forc_dates)
    if (!is.null(object$spec$transform$lambda)) {
        for (i in 1:dim(Y)[3]) {
            Y[,,i] <- box_cox_inverse(Y[,,i], object$spec$transform$lambda[i])
        }
    }
    original_errors <- residuals(object, raw = TRUE)
    out <- lapply(1:NCOL(object$fitted), function(i){
       p <- Y[,,i]
       if (NCOL(p) == 1) p <- matrix(p, ncol = 1)
       colnames(p) <- forc_dates
       class(p) <- "tsmodel.distribution"
       attr(p, "date_class") <- date_class
       p <- list(original_series = zoo(object$spec$target$y_orig[,i], object$spec$target$index), distribution = p)
       class(p) <- "tsmodel.predict"
       error_i <- E[,,i]
       colnames(error_i) <- forc_dates
       class(error_i) <- "tsmodel.distribution"
       L <- list(original_series = as.zoo(original_errors[,i]), distribution = error_i)
       class(L) <- "tsmodel.predict"
       td <- .tsdecompose_predict(States, i, object, forc_dates, newxreg = X[,-1])
       data.table(series = object$spec$target$y_names[i], Level = list(td$Level), 
                  Slope = list(td$Slope), Seasonal = list(td$Seasonal), 
                  X = list(td$X), Error = list(L), Predicted = list(p))
    })
    out <- rbindlist(out)
    L <- list(prediction_table = out, nsim = nsim, h = h, spec = object$spec)
    class(L) <- "tsvets.predict"
    return(L)
}
