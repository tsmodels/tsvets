tsprofile.tsvets.estimate <- function(object, h = 1, nsim = 100, seed = NULL, trace = FALSE, solver = "nlminb", ...)
{
    sim <- simulate(object, seed = seed, nsim = nsim, h = NROW(object$spec$target$y_orig) + h)
    profile <- profile_fun(sim, object, h, trace = trace, solver = solver)
    return(profile)
}

profile_fun <- function(sim, object, h, cores, trace, solver)
{
    n <- NCOL(object$spec$target$y)
    if (trace) {
        prog_trace <- progressor(nrow(sim$simulation_table))
    }
    y_names <- object$spec$target$y_names
    prof %<-% future_lapply(1:nrow(sim$simulation_table), function(i){
        if (trace) prog_trace()
        sim_y <- do.call(cbind, lapply(1:n, function(j) sim$simulation_table$Simulated[[j]]$distribution[j,]))
        colnames(sim_y) <- sim$simulation_table$series
        y <- xts(sim_y, as.POSIXct(rownames(sim_y)))
        yin <- y[1:(nrow(y) - h)]
        spec <- tsspec(object, yin)
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
    }, future.packages = c("tsmethods","tsvets","xts","data.table"), future.seed = TRUE)
    prof <- eval(prof)
    C <- rbindlist(lapply(1:length(prof), function(i) prof[[i]]$L1))
    M <- rbindlist(lapply(1:length(prof), function(i) prof[[i]]$L2))
    
    Actual <- NULL
    Predicted <- NULL
    Simulation <- NULL
    # create distribution for all performance metrics
    stats_distribution <- M[,list(MAPE = mape(Actual, Predicted), MSLRE = mslre(Actual, Predicted), 
                                 BIAS = bias(Actual, Predicted)), by = c("Series","Horizon","Simulation")]
    L <- list(coef = C, true.coef = coef(object), stats_table = stats_distribution)
    return(L)
}
