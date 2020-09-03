########################################################
# A matrix
ss_mat_A <- function(object, type, counter = 0)
{
    switch(type,
           "constant" = ss_mat_A_constant(object, counter),
           "diagonal" = ss_mat_A_diagonal(object, counter),
           "common" = ss_mat_A_common(object, counter),
           "full" = ss_mat_A_full(object, counter),
           "grouped" = ss_mat_A_grouped(object, counter))
    
}

ss_mat_A_common <- function(object, counter = 0)
{
    out <- diag(as.numeric(NA), object$n, object$n)
    N <- object$n * object$n
    parnames <- as.numeric(out)
    parnames[which(is.na(parnames))] <- paste0("Level[Common]")
    parnames[which(!is.na(as.numeric(out)))] <- NA
    parindex <- rep(0, N)
    parindex[which(is.na(as.numeric(out)))] <- counter + 1L
    estimate = rep(0, N)
    estimate[1] <- 1
    lower <- rep(0, N)
    lower[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 1e-12
    upper <- rep(0, N)
    upper[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 1 - 1e-12
    data.table(matrix = "A", parameters = as.numeric(out), parnames = parnames,  index = parindex, estimate = estimate, lower = lower, upper = upper)
}

ss_mat_A_constant <- function(object, counter = 0)
{
    out <- matrix(0, object$n, object$n)
    N <- object$n * object$n
    estimate = rep(0, N)
    lower <- rep(0, N)
    upper <- rep(0, N)
    data.table(matrix = "A", parameters = as.numeric(out), parnames = NA,  index = NA, estimate = estimate, 
               lower = lower, upper = upper)
}

ss_mat_A_diagonal <- function(object, counter = 0)
{
    out <- diag(as.numeric(NA), object$n, object$n)
    N <- object$n * object$n
    parnames <- as.numeric(out)
    parnames[which(is.na(parnames))] <- paste0("Level[",object$series_names,"]")
    parnames[which(!is.na(as.numeric(out)))] <- NA
    parindex <- rep(0, N)
    parindex[which(is.na(as.numeric(out)))] <- counter + (1L:object$n)
    estimate <- rep(0, N)
    estimate[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 1L
    lower <- rep(0, N)
    lower[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 1e-12
    upper <- rep(0, N)
    upper[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 1 - 1e-12
    data.table(matrix = "A", parameters = as.numeric(out), parnames = parnames,  index = parindex, estimate = estimate, 
               lower = lower, upper = upper)
}

ss_mat_A_full <- function(object, counter = 0)
{
    out <- matrix(as.numeric(NA), object$n, object$n)
    N <- object$n * object$n
    parnames <- as.numeric(out)
    parnames[which(is.na(parnames))] <- vectorize_matrix_names(object$series_names, object$series_names)
    parindex <- counter + (1L:N)
    estimate <- rep(1L, N)
    lower <- rep(-1, N)
    lower[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 1e-12
    upper <- rep(1 - 1e-12, N)
    data.table(matrix = "A", parameters = as.numeric(out), parnames = parnames,  index = parindex, estimate = estimate, 
               lower = lower, upper = upper)
}

ss_mat_A_grouped <- function(object, counter = 0)
{
    out <- diag(as.numeric(NA), object$n, object$n)
    N <- object$n * object$n
    parnames <- as.numeric(out)
    groups <- object$group
    ngroups <- length(unique(groups))
    gnames <- paste0("Level[Group = ",groups,"]")
    parnames[which(is.na(parnames))] <- gnames
    parindex <- rep(0, N)
    parindex[which(is.na(as.numeric(out)))] <- counter + groups
    estimate <- rep(0, N)
    ix <- 1 + 0L:(object$n - 1L) * (object$n + 1L)    
    estimate[ix[match(unique(groups), groups)]] <- 1L
    lower <- rep(0, N)
    lower[ix] <- 1e-12
    upper <- rep(0, N)
    upper[ix] <- 1 - 1e-12
    data.table(matrix = "A", parameters = as.numeric(out), parnames = parnames,  index = parindex, estimate = estimate, 
               lower = lower, upper = upper)
}
########################################################
# B matrix

ss_mat_B <- function(object, type, counter = 0)
{
    switch(type,
           "constant" = ss_mat_B_constant(object, counter),
           "diagonal" = ss_mat_B_diagonal(object, counter),
           "common" = ss_mat_B_common(object, counter),
           "full" = ss_mat_B_full(object, counter),
           "grouped" = ss_mat_B_grouped(object, counter))
    
}


ss_mat_B_constant <- function(object, counter = 0)
{
    out <- matrix(0, object$n, object$n)
    N <- object$n * object$n
    estimate = rep(0, N)
    lower <- rep(0, N)
    upper <- rep(0, N)
    data.table(matrix = "B", parameters = as.numeric(out), parnames = NA,  index = NA, estimate = estimate, 
               lower = lower, upper = upper)
}

ss_mat_B_common <- function(object, counter = 0)
{
    out <- diag(as.numeric(NA), object$n, object$n)
    N <- object$n * object$n
    parnames <- as.numeric(out)
    parnames[which(is.na(parnames))] <- paste0("Slope[Common]")
    parnames[which(!is.na(as.numeric(out)))] <- NA
    parindex <- rep(0, N)
    parindex[which(is.na(as.numeric(out)))] <- counter + 1L
    estimate = rep(0, N)
    estimate[1] <- 1
    lower <- rep(0, N)
    lower[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 1e-12
    upper <- rep(0, N)
    upper[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 1 - 1e-12
    data.table(matrix = "B", parameters = as.numeric(out), parnames = parnames,  index = parindex, estimate = estimate, 
               lower = lower, upper = upper)
}

ss_mat_B_diagonal <- function(object, counter = 0)
{
    out <- diag(as.numeric(NA), object$n, object$n)
    N <- object$n * object$n
    parnames <- as.numeric(out)
    parnames[which(is.na(parnames))] <- paste0("Slope[",object$series_names,"]")
    parnames[which(!is.na(as.numeric(out)))] <- NA
    parindex <- rep(0, N)
    parindex[which(is.na(as.numeric(out)))] <- counter + (1L:object$n)
    estimate <- rep(0, N)
    estimate[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 1L
    lower <- rep(0, N)
    lower[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 1e-12
    upper <- rep(0, N)
    upper[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 1 - 1e-12
    data.table(matrix = "B", parameters = as.numeric(out), parnames = parnames,  index = parindex, estimate = estimate, 
               lower = lower, upper = upper)
}


ss_mat_B_full <- function(object, counter = 0)
{
    out <- matrix(as.numeric(NA), object$n, object$n)
    N <- object$n * object$n
    parnames <- as.numeric(out)
    parnames[which(is.na(parnames))] <- vectorize_matrix_names(object$series_names, object$series_names, state_name = "Slope")
    parindex <- counter + (1L:N)
    estimate <- rep(1L, N)
    lower <- rep(-1, N)
    lower[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 1e-12
    upper <- rep(1 - 1e-12, N)
    data.table(matrix = "B", parameters = as.numeric(out), parnames = parnames,  index = parindex, estimate = estimate, 
               lower = lower, upper = upper)
}


ss_mat_B_grouped <- function(object, counter = 0)
{
    out <- diag(as.numeric(NA), object$n, object$n)
    N <- object$n * object$n
    parnames <- as.numeric(out)
    groups <- object$group
    ngroups <- length(unique(groups))
    gnames <- paste0("Slope[Group = ",groups,"]")
    parnames[which(is.na(parnames))] <- gnames
    parindex <- rep(0, N)
    parindex[which(is.na(as.numeric(out)))] <- counter + groups
    estimate <- rep(0, N)
    ix <- 1 + 0L:(object$n - 1L) * (object$n + 1L)    
    estimate[ix[match(unique(groups), groups)]] <- 1L
    lower <- rep(0, N)
    lower[ix] <- 1e-12
    upper <- rep(0, N)
    upper[ix] <- 1 - 1e-12
    data.table(matrix = "B", parameters = as.numeric(out), parnames = parnames,  index = parindex, estimate = estimate, 
               lower = lower, upper = upper)
}
########################################################
# Phi Matrix

ss_mat_Phi <- function(object, type = "none", counter)
{
    switch(type,
           "none" = ss_mat_Phi_default(object, counter),
           "diagonal" = ss_mat_Phi_diagonal(object, counter),
           "common" = ss_mat_Phi_common(object, counter),
           "full" = ss_mat_Phi_full(object, counter),
           "grouped" = ss_mat_Phi_grouped(object, counter))
}

ss_mat_Phi_default <- function(object, counter = 0)
{
    out <- diag(1, object$n, object$n)
    parnames <- as.numeric(out)
    #data.table(matrix = "Phi", parameters = as.numeric(out), parnames = NA,  index = NA, estimate = 0, lower = 1, upper = 1)
    data.table(matrix = "Phi", parameters = as.numeric(out), parnames = NA,  index = 0, estimate = 0, lower = 1, upper = 1)
}

ss_mat_Phi_common <- function(object, counter = 0)
{
    out <- diag(as.numeric(NA), object$n, object$n)
    N <- object$n * object$n
    parnames <- as.numeric(out)
    parnames[which(is.na(parnames))] <- paste0("Dampening[Common]")
    parnames[which(!is.na(as.numeric(out)))] <- NA
    parindex <- rep(0, N)
    parindex[which(is.na(as.numeric(out)))] <- counter + 1L
    estimate = rep(0, N)
    estimate[1] <- 1
    lower <- rep(0, N)
    lower[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 0.6
    upper <- rep(0, N)
    upper[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 1
    data.table(matrix = "Phi", parameters = as.numeric(out), parnames = parnames,  index = parindex, estimate = estimate, 
               lower = lower, upper = upper)
}

ss_mat_Phi_diagonal <- function(object, counter = 0)
{
    out <- diag(as.numeric(NA), object$n, object$n)
    N <- object$n * object$n
    parnames <- as.numeric(out)
    parnames[which(is.na(parnames))] <- paste0("Dampening[",object$series_names,"]")
    parnames[which(!is.na(as.numeric(out)))] <- NA
    parindex <- rep(0, N)
    parindex[which(is.na(as.numeric(out)))] <- counter + (1L:object$n)
    estimate <- rep(0, N)
    estimate[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 1L
    lower <- rep(0, N)
    lower[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 0.6
    upper <- rep(0, N)
    upper[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 1
    data.table(matrix = "Phi", parameters = as.numeric(out), parnames = parnames,  index = parindex, estimate = estimate, 
               lower = lower, upper = upper)
}


ss_mat_Phi_full <- function(object, counter = 0)
{
    out <- matrix(as.numeric(NA), object$n, object$n)
    N <- object$n * object$n
    parnames <- as.numeric(out)
    parnames[which(is.na(parnames))] <- vectorize_matrix_names(object$series_names, object$series_names, state_name = "Dampening")
    parindex <- counter + (1L:N)
    estimate <- rep(1L, N)
    lower <- rep(-1, N)
    lower[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 0.6
    upper <- rep(1 - 1e-12, N)
    data.table(matrix = "Phi", parameters = as.numeric(out), parnames = parnames,  index = parindex, estimate = estimate, 
               lower = lower, upper = upper)
}


ss_mat_Phi_grouped <- function(object, counter = 0)
{
    out <- diag(as.numeric(NA), object$n, object$n)
    N <- object$n * object$n
    parnames <- as.numeric(out)
    groups <- object$group
    ngroups <- length(unique(groups))
    gnames <- paste0("Dampening[Group = ",groups,"]")
    parnames[which(is.na(parnames))] <- gnames
    parindex <- rep(0, N)
    parindex[which(is.na(as.numeric(out)))] <- counter + groups
    estimate <- rep(0, N)
    ix <- 1 + 0L:(object$n - 1L) * (object$n + 1L)    
    estimate[ix[match(unique(groups), groups)]] <- 1L
    lower <- rep(0, N)
    lower[ix] <- 1e-12
    upper <- rep(0, N)
    upper[ix] <- 1 - 1e-12
    data.table(matrix = "Phi", parameters = as.numeric(out), parnames = parnames,  index = parindex, estimate = estimate, 
               lower = lower, upper = upper)
}

########################################################
# K Matrix
ss_mat_K <- function(object, type, counter = 0)
{
    switch(type,
           "diagonal" = ss_mat_K_diagonal(object, counter),
           "common" = ss_mat_K_common(object, counter),
           "full" = ss_mat_K_full(object, counter),
           "grouped" = ss_mat_K_grouped(object, counter))
    
}

ss_mat_K_common <- function(object, counter = 0)
{
    out <- diag(as.numeric(NA), object$n, object$n)
    N <- object$n * object$n
    parnames <- as.numeric(out)
    parnames[which(is.na(parnames))] <- paste0("Seasonal[Common]")
    parnames[which(!is.na(as.numeric(out)))] <- NA
    parindex <- rep(0, N)
    parindex[which(is.na(as.numeric(out)))] <- counter + 1L
    estimate = rep(0, N)
    estimate[1] <- 1
    lower <- rep(0, N)
    lower[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 1e-12
    upper <- rep(0, N)
    upper[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 1 - 1e-12
    data.table(matrix = "K", parameters = as.numeric(out), parnames = parnames,  index = parindex, estimate = estimate, 
               lower = lower, upper = upper)
}

ss_mat_K_diagonal <- function(object, counter = 0)
{
    out <- diag(as.numeric(NA), object$n, object$n)
    N <- object$n * object$n
    parnames <- as.numeric(out)
    parnames[which(is.na(parnames))] <- paste0("Seasonal[",object$series_names,"]")
    parnames[which(!is.na(as.numeric(out)))] <- NA
    parindex <- rep(0, N)
    parindex[which(is.na(as.numeric(out)))] <- counter + (1L:object$n)
    estimate <- rep(0, N)
    estimate[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 1L
    lower <- rep(0, N)
    lower[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 1e-12
    upper <- rep(0, N)
    upper[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 1 - 1e-12
    data.table(matrix = "K", parameters = as.numeric(out), parnames = parnames,  index = parindex, estimate = estimate, 
               lower = lower, upper = upper)
}


ss_mat_K_full <- function(object, counter = 0)
{
    out <- matrix(as.numeric(NA), object$n, object$n)
    N <- object$n * object$n
    parnames <- as.numeric(out)
    parnames[which(is.na(parnames))] <- vectorize_matrix_names(object$series_names, object$series_names, state_name = "Seasonal")
    parindex <- counter + (1L:N)
    estimate <- rep(1L, N)
    lower <- rep(-1, N)
    lower[1 + 0L:(object$n - 1L) * (object$n + 1L)] <- 1e-12
    upper <- rep(1 - 1e-12, N)
    data.table(matrix = "K", parameters = as.numeric(out), parnames = parnames,  index = parindex, estimate = estimate, 
               lower = lower, upper = upper)
}


ss_mat_K_grouped <- function(object, counter = 0)
{
    out <- diag(as.numeric(NA), object$n, object$n)
    N <- object$n * object$n
    parnames <- as.numeric(out)
    groups <- object$group
    ngroups <- length(unique(groups))
    gnames <- paste0("Seasonal[Group = ",groups,"]")
    parnames[which(is.na(parnames))] <- gnames
    parindex <- rep(0, N)
    parindex[which(is.na(as.numeric(out)))] <- counter + groups
    estimate <- rep(0, N)
    ix <- 1 + 0L:(object$n - 1L) * (object$n + 1L)    
    estimate[ix[match(unique(groups), groups)]] <- 1L
    lower <- rep(0, N)
    lower[ix] <- 1e-12
    upper <- rep(0, N)
    upper[ix] <- 1 - 1e-12
    data.table(matrix = "K", parameters = as.numeric(out), parnames = parnames,  index = parindex, estimate = estimate, 
               lower = lower, upper = upper)
}
########################################################
# xreg
ss_mat_xreg <- function(object, counter)
{
    # xin <- matrix(c(0,1,2,2,1,0,1,2,1,1,2,2), ncol = 4, nrow = 3, byrow = TRUE)
    n <- object$n
    m <- ncol(object$xreg)
    xin <- object$xreg_include
    if (ncol(xin) != m) stop("\nncol xreg_include not equal to ncol xreg")
    if (nrow(xin) != n) stop("\nnrow xreg_include not equal to ncol y")
    if (all(xin == 0)) stop("\nall values of xreg_include are zero")
    k <- counter
    env <- environment()
    env$k <- k
    xp <- lapply(1:ncol(xin), function(i){
        k <- env$k
        estimate <- rep(0, n)
        index <- rep(0, n)
        xmat <- matrix(NA, ncol = 1, nrow = n)
        use <- which(xin[,i] == 1)
        if (length(use) > 0 ) {
            xmat[use, 1] <- paste0("Beta[",use,"_",i,"]")
            estimate[use] <- 1
            index[use] <- k + 1:length(use)
            assign("k",k + length(use),envir = env)
        }
        k <- env$k
        use <- which(xin[,i] == 2)
        if (length(use) > 0 ) {
            xmat[use, 1] <- paste0("Beta[",i,"]")
            estimate[use[1]] <- 1
            index[use] <- k + 1
            assign("k",k + 1,envir = env)
        }
        return(list(xmat = xmat, estimate = estimate, index = index))
    })
    xnames <- do.call(cbind, lapply(1:length(xp), function(i) xp[[i]]$xmat))
    estimate <- do.call(c, lapply(1:length(xp), FUN = function(i) c(xp[[i]]$estimate)))
    index <- do.call(c, lapply(1:length(xp), FUN = function(i) c(xp[[i]]$index)))
    pars <- as.numeric(xin)
    pars[which(pars > 0)] <- NA
    lower <- rep(0, length(pars))
    lower[which(is.na(pars))] <- -100
    upper <- rep(0, length(pars))
    upper[which(is.na(pars))] <-  100
    data.table(matrix = "X", parameters = pars, parnames = as.vector(xnames), index = index, estimate = estimate, 
               lower = lower, upper = upper)
}
ss_mat_H <- function(object)
{
    H <- NULL
    D <- Matrix(diag(1, object$n, object$n), sparse = TRUE)
    H <- rbind(H, D) 
    if (object$slope) {
        H <- rbind(H, D)
    }
    if (object$seasonal) {
        H <- rbind(H, Matrix(0, nrow = object$n * (object$frequency - 1), ncol = object$n, sparse = TRUE))
        H <- rbind(H, D)
    }
    return(H)
}
ss_mat_Fbar <- function(object)
{
    cbind(rbind(Matrix(rep(0, object$frequency - 1), nrow = 1, sparse = TRUE),
          Matrix(diag(1, object$frequency - 1, object$frequency - 1), sparse = TRUE)), 
          Matrix(c(1, rep(0, object$frequency - 1)), nrow = object$frequency, sparse = TRUE))
}

ss_mat_F <- function(object, Phi)
{
    if (object$slope & object$seasonal) {
        D <- Matrix(diag(object$n), sparse = TRUE)
        K <-  kronecker(ss_mat_Fbar(object), D)
        O <- Matrix(0, object$n, object$n, sparse = TRUE)
        M <- Matrix(0, object$n * object$frequency, object$n, sparse = TRUE)
        FF <- rbind(cbind(D, D, t(M)), cbind(O, Phi, t(M)), cbind(M, M, K))
    } else if (!object$slope & object$seasonal) {
        D <- Matrix(diag(object$n), sparse = TRUE)
        K <-  kronecker(ss_mat_Fbar(object), D)
        M <- Matrix(0, object$n * object$frequency, object$n, sparse = TRUE)
        FF <- rbind(cbind(D, t(M)), cbind(M, K))
    } else if (object$slope & !object$seasonal) {
        D <- Matrix(diag(object$n), sparse = TRUE)
        O <- Matrix(0, object$n, object$n, sparse = TRUE)
        M <- Matrix(0, object$n * object$frequency, object$n, sparse = TRUE)
        FF <- rbind(cbind(D, D), cbind(O, Phi))
    } else if (!object$slope & !object$seasonal) {
        FF <- Matrix(diag(object$n), sparse = TRUE)
    } else {
        stop("\nunrecognized combination of states (Fmat)")
    }
    return(FF)
}
    

ss_mat_Gbar <- function(object)
{
    Matrix(diag(c( (object$frequency - 1)/object$frequency, rep(-1/object$frequency, object$frequency - 1)), 
                object$frequency, object$frequency), sparse = TRUE)
}

ss_mat_G <- function(object)
{
    if (object$seasonal) {
        if (object$slope) {
            D <- Matrix(diag(object$n), sparse = TRUE)
            K <-  kronecker(ss_mat_Gbar(object), D)
            O <- Matrix(0, nrow = object$n, ncol = object$n, sparse = TRUE)
            M <- Matrix(0, nrow = object$n * object$frequency, ncol = object$n, sparse = TRUE)
            G <- rbind(cbind(D, O, t(M)), cbind(O, D, t(M)), cbind(M, M, K))
        } else {
            D <- Matrix(diag(object$n), sparse = TRUE)
            K <-  kronecker(ss_mat_Gbar(object), D)
            M <- Matrix(0, nrow = object$n * object$frequency, ncol = object$n, sparse = TRUE)
            G <- rbind(cbind(D, t(M)), cbind(M, K))
        }
    } else {
        if (object$slope) {
            D <- Matrix(diag(object$n), sparse = TRUE)
            O <- Matrix(0, nrow = object$n, ncol = object$n, sparse = TRUE)
            G <- rbind(cbind(D, O), cbind(O, D))
        } else {
            G <- Matrix(diag(object$n), sparse = TRUE)
        }
    }
    return(G)
}

vectorize_matrix_names <- function(row_names, col_names, state_name = "Level")
{
    # need to reverse switch the columns to give column wise vectorization
    g <- expand.grid(col_names, row_names)[,c(2,1)]
    paste0(state_name,"[",g[,1], "_",g[,2],"]")
}

bdiag <- function(x = NULL, ...)
{
    if (is.null(x)) {
        if (nargs() == 1L) x <- as.list(...) else x <- list(...)
    }
    n <- length(x)
    if (n == 0L) return(NULL)
    x <- lapply(x, function(y) if (length(y)) as.matrix(y) else stop("Zero-length component in x"))
    d <- array(unlist(lapply(x, dim)), c(2, n))
    rr <- d[1L, ]
    cc <- d[2L, ]
    rsum <- sum(rr)
    csum <- sum(cc)
    out <- array(0, c(rsum, csum))
    ind <- array(0, c(4, n))
    rcum <- cumsum(rr)
    ccum <- cumsum(cc)
    ind[1L, -1L] <- rcum[-n]
    ind[2L, ] <- rcum
    ind[3L, -1L] <- ccum[-n]
    ind[4L, ] <- ccum
    imat <- array(1:(rsum * csum), c(rsum, csum))
    iuse <- apply(ind, 2, function(y, imat) imat[(y[1L] + 1):y[2L], (y[3L] + 1):y[4L]], imat = imat)
    iuse <- as.vector(unlist(iuse))
    out[iuse] <- unlist(x)
    return(out)
}

p_matrix <- function(object)
{
    Amat <- object$spec$vets_env$Amat
    n <- NCOL(object$spec$target$y_orig)
    Level_matrix <- Matrix(t(Amat[, 1:n, drop = FALSE]))
    k <- n
    if (object$spec$model$slope != "none") {
        Slope_matrix <- Matrix(t(Amat[, (k + 1):(k + n), drop = FALSE]))
        if (object$spec$model$damped != "none") {
            Phi_matrix <- Matrix(t(Amat[, object$spec$vets_env$Phi_index[1]:object$spec$vets_env$Phi_index[2], drop = FALSE]))
        } else {
            Phi_matrix <- NULL
        }
        k <- k + n
    } else {
        Slope_matrix <- NULL
        Phi_matrix <- NULL
    }
    if (object$spec$model$seasonal != "none") {
        Seasonal_matrix <- Matrix(t(Amat[, (k + 1):(k + n)]))
    } else {
        Seasonal_matrix <- NULL
    }
    if (object$spec$xreg$include_xreg) {
        X_matrix <- Matrix(Amat[,object$spec$vets_env$X_index[1]:object$spec$vets_env$X_index[2], drop = FALSE])
        colnames(X_matrix) <- colnames(object$spec$xreg$xreg)
        rownames(X_matrix) <- object$spec$target$y_names
    } else {
        X_matrix <- NULL
    }
    return(list(Level_matrix = Level_matrix, Slope_matrix = Slope_matrix, Phi_matrix = Phi_matrix, Seasonal_matrix = Seasonal_matrix, X_matrix = X_matrix))
}