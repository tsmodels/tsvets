tsprofile.tsvets.estimate <- function(object, h = 1, nsim = 100, seed = NULL, cores = 1, trace = 0, solver = "nlminb", ...)
{
    sim <- simulate(object, seed = seed, nsim = nsim, h = NROW(object$spec$target$y_orig) + h)
    profile <- profile_fun(sim, object, h, cores = cores, trace = trace, solver = solver)
    return(profile)
}

profile_fun <- function(sim, object, h, cores, trace, solver)
{
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    if (trace == 1) {
        iterations <- nrow(sim$simulation_table)
        pb <- txtProgressBar(max = iterations, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
    } else {
        opts <- NULL
    }
    n <- NCOL(object$spec$target$y)
    i <- 1
    y_names <- object$spec$target$y_names
    P <- foreach(i = 1:nrow(sim$simulation_table), .packages = c("tsmethods","tsvets","xts","data.table"), .options.snow = opts) %dopar% {
        sim_y <- do.call(cbind, lapply(1:n, function(j) sim$simulation_table$Simulated[[j]]$distribution[j,]))
        colnames(sim_y) <- sim$simulation_table$series
        y <- xts(sim_y, as.POSIXct(rownames(sim_y)))
        yin <- y[1:(nrow(y) - h)]
        spec <- tsspec(object, yin, lambda = object$spec$transform$lambda)
        # add try catch
        mod <- try(estimate(spec, solver = solver), silent = TRUE)
        if (inherits(mod, 'try-error')) {
            return(list(L1 = NULL, L2 = NULL))
        }
        p <- predict(mod, h = h)
        L1 <- data.table("Variable" = names(coef(mod)), "Value" = coef(mod), "Simulation" = i)
        L2 <- lapply(1:n, function(j){
            data.table(Series = y_names[j], "Predicted" = as.numeric(colMeans(p$prediction_table[j]$Predicted[[1]]$distribution)), 
                       "Actual" = as.numeric(tail(y[,j], h)), "Simulation" = i, "Horizon" = 1:h)
        })
        L2 <- rbindlist(L2)
        return(list(L1 = L1, L2 = L2))
    }
    C <- rbindlist(lapply(1:length(P), function(i) P[[i]]$L1))
    M <- rbindlist(lapply(1:length(P), function(i) P[[i]]$L2))
    
    if (trace == 1) {
        close(pb)
    }
    stopCluster(cl)
    Actual <- NULL
    Predicted <- NULL
    Simulation <- NULL
    # create distribution for all performance metrics
    stats_distribution <- M[,list(MAPE = mape(Actual, Predicted), MSLRE = mslre(Actual, Predicted), 
                                 BIAS = bias(Actual, Predicted)), by = c("Series","Horizon","Simulation")]
    coef_distribution <- C
    
    L <- list(coef = C, true.coef = coef(object), stats_table = stats_distribution)
    return(L)
}
