#' Series Aggregation
#'
#' @description Aggregates estimated, predicted or simulated series based on a 
#' weighting vector.
#' @param object object of class \dQuote{tsvets.estimate}. \dQuote{tsvets.predict} 
#' or \dQuote{tsvets.simulate}.
#' @param weights vector of weights of length equal to the number of series.
#' @param return_model if the estimated object is a homogeneous coefficients model 
#' with common lamda parameter (of either 0, 1 or NULL), then it will return a 
#' univariate tsvets object of class \dQuote{tsvets.estimate}.
#' @param ... not currently used.
#' @return Depends on the input class.
#' @details For an estimated object which has common components (homogeneous 
#' coefficients) for all estimated states and common lambda parameter, then a 
#' reconstructed object of the same class representing the weighted 
#' representation of the model is returned. In all other cases and input classes, 
#' the returned object will depend on whether the lambda parameter was 0 or 1 
#' for all underlying series. For the case of the log transform (lambda = 0), 
#' then the states aggregate (given a weighting vector) whilst the actual, fitted, 
#' predicted or simulated values aggregate in logs and are then exponentiated 
#' (a multiplicative model). In all other cases only the actual, fitted, predicted 
#' or simulated values are returned representing the weighted aggregation of the 
#' underlying series given the weighting vector.
#' @references Hyndman, Rob and Koehler, Anne B and Ord, J Keith and Snyder, 
#' Ralph D, 2008, Forecasting with exponential smoothing: the state space approach, 
#' Section 17.1.2,,Springer Science \& Business Media.
#' @aliases tsaggregate
#' @method tsaggregate tsvets.estimate
#' @rdname tsaggregate
#' @export
#'
#'
#'
tsaggregate.tsvets.estimate <- function(object, weights = NULL, return_model = FALSE, ...)
{
    if (is.null(object$spec$transform[[1]])) {
        lambda <- rep(1, ncol(object$fitted))
    } else {
        lambda <- sapply(object$spec$transform, function(x) x$lambda)
    }
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
                               seasonal = object$spec$model$seasonal, frequency = object$spec$target$frequency, 
                               dependence = "diagonal", lambda = lambda[1], 
                               xreg = xreg, xreg_include = xreg_include)
        xseed <- object$spec$vets_env$init_states %*% weights
        n_states <- NROW(spec$vets_env$States)
        # adjust the initial states in the Amat matrix
        spec$vets_env$Amat[,(ncol(spec$vets_env$Amat) - n_states + 1):ncol(spec$vets_env$Amat)] <- xseed[1:n_states]
        spec$vets_env$States[,1] <- xseed[1:NROW(spec$vets_env$States)]
        spec$vets_env$good <- t(spec$target$good_matrix)
        spec$vets_env$selection <- spec$target$good_index
        spec$vets_env$ymat <- na.fill(spec$vets_env$ymat, fill =  0)
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
        states <- do.call(rbind, lapply(1:nrow(object$States), function(i){
            matrix(matrix(object$States[i,], ncol = n, byrow = T) %*% weights, nrow = 1)
        }))
        if (!is.null(spec$transform)) {
            fit <- do.call(cbind, lapply(1:ncol(object$fitted), function(i) spec$transform[[1]]$transform(object$fitted[,i], spec$transform[[1]]$lambda)))
            fit <- fit %*% weights
            e <- spec$transform[[1]]$transform(actual,spec$transform[[1]]$lambda)  - fit
            fit <- xts(spec$transform[[1]]$inverse(fit, spec$transform[[1]]$lambda), spec$target$index)
        } else {
            fit <- object$fitted %*% weights
            e <- actual - fit
        }
        spec$vets_env$Amat <- filt$Amat
        spec$vets_env$Fmat <- filt$Fmat
        spec$vets_env$Gmat <- filt$Gmat
        out <- list(fitted = fit, States = states, Error = e, spec = spec, opt = list(pars = pars, negative_llh = NA), variance = V)
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
                if (all(lambda == 0)) {
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
                if (!is.null(object$spec$transform)) {
                    if (all(lambda == 0)) {
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

#' @method tsaggregate tsvets.predict
#' @rdname tsaggregate
#' @export
tsaggregate.tsvets.predict <- function(object, weights = NULL, ...)
{
    n <- NCOL(object$spec$target$y)
    if (is.null(weights)) {
        weights <- matrix(1, ncol = n, nrow = 1)
    } else {
        if (length(as.numeric(weights)) != n) stop("\nweights should be a vector of length equal to ncol y.")   
        weights <- matrix(as.numeric(weights), ncol = 1, nrow = n)
    }
    if (is.null(object$spec$transform[[1]])) {
        lambda <- rep(1, NROW(object$prediction_table))
    } else {
        lambda <- sapply(object$spec$transform, function(x) x$lambda)
    }
    condition_transform <- (all(lambda == 0) | all(lambda == 1))
    condition_homogeneous <- object$spec$model$level == "common" & object$spec$model$slope %in% c("none","common") & object$spec$model$seasonal %in% c("none","common") & object$spec$model$damped %in% c("none","common")
    
    if (!condition_homogeneous) {
        if (condition_transform) {
            f_dates <- colnames(object$prediction_table[1]$Predicted[[1]]$distribution)
            date_class <- attr(object$prediction_table[1]$Predicted[[1]]$distribution, "date_class")
            array_predicted <- array(unlist(lapply(1:nrow(object$prediction_table), function(i) as.matrix(object$prediction_table[i]$Predicted[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$prediction_table)))
            array_predicted <- aperm(array_predicted, perm = c(1,3,2))
            if (all(lambda == 0)) {
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
            if (!is.null(object$spec$transform[[1]])) {
                # we can aggregate states for lambda = 0, 1, or NULL
                if (all(lambda == 0)) {
                    f_dates <- colnames(object$prediction_table[1]$Predicted[[1]]$distribution)
                    date_class <- attr(object$prediction_table[1]$Predicted[[1]]$distribution, "date_class")
                    array_predicted <- array(unlist(lapply(1:nrow(object$prediction_table), function(i) as.matrix(object$prediction_table[i]$Predicted[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$prediction_table)))
                    array_predicted <- aperm(array_predicted, perm = c(1,3,2))
                    distribution <- apply(array_predicted, 3, function(x) exp(log(x) %*% weights))
                    colnames(distribution) <- f_dates
                    class(distribution) <- "tsmodel.distribution"
                    attr(distribution, "date_class") <- date_class
                    original_series <- xts(exp(log(object$spec$target$y_orig) %*% weights), object$spec$target$index)
                } else if (all(lambda == 1)) {
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

#' @method tsaggregate tsvets.simulate
#' @rdname tsaggregate
#' @export
tsaggregate.tsvets.simulate <- function(object, weights = NULL, ...)
{
    n <- NCOL(object$spec$target$y)
    if (is.null(weights)) {
        weights <- matrix(1, ncol = n, nrow = 1)
    } else {
        if (length(as.numeric(weights)) != n) stop("\nweights should be a vector of length equal to ncol y.")   
        weights <- matrix(as.numeric(weights), ncol = 1, nrow = n)
    }
    if (is.null(object$spec$transform[[1]])) {
        lambda <- rep(1, NROW(object$simulation_table))
    } else {
        lambda <- sapply(object$spec$transform, function(x) x$lambda)
    }
    condition_transform <- (all(lambda == 0) | all(lambda == 1))
    condition_homogeneous <- object$spec$model$level == "common" & object$spec$model$slope %in% c("none","common") & object$spec$model$seasonal %in% c("none","common") & object$spec$model$damped %in% c("none","common")
    
    if (!condition_homogeneous) {
        if (condition_transform) {
            f_dates <- colnames(object$simulation_table[1]$Simulated[[1]]$distribution)
            date_class <- attr(object$simulation_table[1]$Simulated[[1]]$distribution, "date_class")
            array_predicted <- array(unlist(lapply(1:nrow(object$simulation_table), function(i) as.matrix(object$simulation_table[i]$Simulated[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$simulation_table)))
            array_predicted <- aperm(array_predicted, perm = c(1,3,2))
            if (all(lambda == 0)) {
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
            if (!is.null(object$spec$transform)) {
                # we can aggregate states for lambda = 0, 1, or NULL
                if (all(lambda == 0)) {
                    f_dates <- colnames(object$simulation_table[1]$Simulated[[1]]$distribution)
                    date_class <- attr(object$simulation_table[1]$Simulated[[1]]$distribution, "date_class")
                    array_predicted <- array(unlist(lapply(1:nrow(object$simulation_table), function(i) as.matrix(object$simulation_table[i]$Simulated[[1]]$distribution))), dim = c(object$nsim, object$h, nrow(object$simulation_table)))
                    array_predicted <- aperm(array_predicted, perm = c(1,3,2))
                    distribution <- apply(array_predicted, 3, function(x) exp(log(x) %*% weights))
                    colnames(distribution) <- f_dates
                    class(distribution) <- "tsmodel.distribution"
                    attr(distribution, "date_class") <- date_class
                    original_series <- xts(exp(log(object$spec$target$y_orig) %*% weights), object$spec$target$index)
                } else if (all(lambda == 1)) {
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