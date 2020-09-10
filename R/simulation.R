simulate.tsvets.estimate <- function(object, nsim = 1, seed = NULL, h = NULL, newxreg = NULL, sim_dates = NULL, bootstrap = FALSE,
                                     pars = coef(object), init_states = object$spec$vets_env$States[,1], ...)
{
    if (is.null(seed)) {
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    } else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    if (is.null(h)) h = nrow(object$spec$target$y_orig)
    xseed_length <- length(object$spec$vets_env$States[,1])
    if (is.null(init_states)) {
        xseed <- object$spec$vets_env$States[,1,drop = FALSE]
    } else {
        if (length(as.numeric(init_states)) != xseed_length) stop(paste0("\ninit.states must be of dimensions 1 x ", xseed_length))
        xseed <- init_states
    }
    xseed <- t(xseed)
    
    cpars <- coef(object)
    if (!is.null(pars)) {
        pars <- na.omit(pars[names(cpars)])
        if (length(pars) == 0) {
            stop("\npars names do not match any of the model parameters")
        }
        cpars[names(pars)] <- pars
    }
    
    if (!object$spec$xreg$include_xreg) {
        X <- matrix(0, ncol = 1, nrow = h + 1)
        if (is.null(sim_dates)) {
            sim_dates = future_dates(tail(object$spec$target$index, 1), frequency = object$spec$target$sampling, n = h)
        }
    } else {
        if (!is.null(newxreg)) {
            sim_dates = index(newxreg)
            X <- rbind(matrix(0, ncol = ncol(newxreg), nrow = 1), coredata(newxreg))
        }
        else {
            if (is.null(sim_dates)) {
                sim_dates = future_dates(tail(object$spec$target$index, 1), frequency = object$spec$target$sampling, n = h)
            }
            warning("\nxreg use in estimation but newxreg is NULL...setting to zero")
            X <- matrix(0, ncol = ncol(object$spec$xreg$xreg), nrow = h + 1)
        }
    }
    X <- t(X)
    date_class <- attr(object$spec$target$sampling, "date_class")
    
    n <- ncol(object$Error)
    innov <- object$Error
    innov_index <- 1:nrow(innov)
    E <- array(0, dim = c(n, h + 1, nsim))
    S <- tscov(object)
    if (bootstrap) {
        for (i in 1:nsim) {
            E[,,i] <- t(innov[sample(innov_index, h + 1, replace = TRUE),,drop = FALSE])
        }
    } else {
        for (i in 1:nsim) {
            E[,,i] <- t(mvrnorm(h + 1, mu = rep(0, n), as.matrix(S)))
        }
    }
    ##################################################################
    Amat <- object$spec$vets_env$Amat
    Amat[object$spec$vets_env$matrix_index] <- cpars[object$spec$vets_env$parameter_index]
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
    model <- c(h, nsim, n, ifelse(object$spec$xreg$include_xreg, 1, 0))
    f <- vets_cpp_simulate(model, E, Amat, Fmat, Hmat, Gmat, istate = xseed, X, beta)
    Y <- aperm(f$Y, map = list(3,2,1))
    States <- aperm(f$States, list(3,2,1))
    sim_dates <- as.character(sim_dates)
    if (!is.null(object$spec$transform$lambda)) {
        for (i in 1:dim(Y)[3]) {
            Y[,,i] <- box_cox_inverse(Y[,,i], object$spec$transform$lambda[i])
        }
    }
    out <- lapply(1:n, function(i){
        p <- Y[,,i]
        if (NCOL(p) == 1) p <- matrix(p, ncol = length(sim_dates))
        colnames(p) <- sim_dates
        class(p) <- c("tsmodel.distribution")
        attr(p, date_class) <- date_class
        p <- list(original_series = NULL, distribution = p)
        class(p) <- c("tsmodel.predict")
        td <- .tsdecompose_simulate(States, i, object, sim_dates)
        data.table(series = object$spec$target$y_names[i], Level = list(td$Level), Slope = list(td$Slope), Seasonal = list(td$Seasonal), Simulated = list(p))
    })
    out <- rbindlist(out)
    L <- list(simulation_table = out, nsim = nsim, h = h, spec = object$spec)
    class(L) <- "tsvets.simulate"
    return(L)
}
