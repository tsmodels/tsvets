tsfilter.tsvets.estimate <- function(object, y = NULL, newxreg = NULL, ...)
{
    yold <- xts(object$spec$target$y_orig, object$spec$target$index)
    ynew <- y
    exc <- which(index(ynew) %in% index(yold))
    if (length(exc) == 0) {
        y <- ynew
    } else{
        y <- ynew[-exc]
        if (NROW(y) == 0) {
            warning("\nno new data in y which is not already in object!")
            return(object)
        }
    }
    if (object$spec$xreg$include_xreg) {
        nx <- NCOL(object$spec$xreg$xreg)
        if (!is.null(newxreg)) {
            if (ncol(newxreg) != nx) stop(paste0("\nExpected ", nx, " columns in newxreg but got ", NCOL(newxreg)))
            newcnames <- colnames(newxreg)
            oldcnames <- colnames(object$spec$xreg$xreg)
            if (!is.null(newcnames) & !is.null(oldcnames)) {
                if (!all(sort(oldcnames) %in% sort(newcnames))) {
                    stop("\ncolnames of newxreg do not match those of original xreg")
                }
                newxreg <- newxreg[, oldcnames]
            }
            X <- coredata(newxreg)
            if (length(exc) > 0) {
                X <- X[-exc,,drop = FALSE]
            }
        } else {
            X <- matrix(0, ncol = nx, nrow = nrow(y))
        }
    } else {
        X <- matrix(0, ncol = 1, nrow = NROW(y))
        if (length(exc) > 0) {
            X <- X[-exc,,drop = FALSE]
        }
    }
    newindex <- index(y)
    yneworig <- y
    if (!is.null(object$spec$transform)) {
        y <- object$spec$transform$transform(y, object$spec$transform$lambda)
    }
    xseed <- tail(object$States, 1)
    pars <- object$opt$par
    env <- new.env()
    for (n in ls(object$spec$vets_env, all.names = TRUE)) assign(n, get(n, object$spec$vets_env), env)
    n <- nrow(ynew)
    env$States <- t(rbind(xseed, matrix(0, ncol = ncol(xseed), nrow = n)))
    env$xreg <- t(rbind(matrix(0, ncol = ncol(X), nrow = 1), X))
    env$model[1] <- n + 1
    env$ymat <- t(rbind(matrix(0, ncol = ncol(y), nrow = 1), coredata(y)))
    f <- vets_filter(pars, env)
    # augment the original object and return back
    if (!is.null(object$spec$transform)) {
        y_fit <- xts(object$spec$transform$inverse(f$fitted[-1,], object$spec$transform$lambda), newindex)
    } else {
        y_fit <- xts(f$fitted[-1,], newindex)
    }
    
    object$States <- rbind(object$States, f$States[-1,])
    object$fitted <- rbind(object$fitted, y_fit)
    object$Error <- rbind(object$Error, f$Error[-1,])
    object$spec$target$y <- rbind(object$spec$target$y, coredata(y))
    object$spec$target$y_orig <- rbind(object$spec$target$y_orig, coredata(yneworig))
    object$spec$target$index <- c(object$spec$target$index, newindex)
    if (object$spec$xreg$include_xreg) {
        object$spec$xreg$xreg <- rbind(object$spec$xreg$xreg, X)
        object$spec$vets_env$xreg <- cbind(object$spec$vets_env$xreg, t(X))
    }
    object$spec$vets_env$model[1] <- nrow(object$fitted) + 1
    object$spec$vets_env$ymat <- cbind(object$spec$vets_env$ymat, t(coredata(y)))
    return(object)
}
