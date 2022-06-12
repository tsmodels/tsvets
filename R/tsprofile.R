#' Model Simulation Based Profiling
#'
#' @description Profiling of model dynamics using simulation/estimation/prediction.
#' @details The function profiles an estimated model by simulating and then 
#' estimating multiple paths from the assumed DGP while leaving h values out for 
#' prediction evaluation. Each simulated path is equal to the size of the original 
#' dataset plus h additional values, and initialized with the initial state vector 
#' from the model. A data.table matrix is returned with the distribution of the 
#' coefficients from each path estimation as well as a stats table with the MAPE, 
#' BIAS and MSLRE by horizon, simulation and series.
#' @param object an object of class \dQuote{tsvets.estimate}.
#' @param h the forecast horizon on which to evaluate performance metrics.
#' @param nsim the number of paths to generate.
#' @param seed an object specifying if and how the random number generator
#' should be initialized. See the simulate documentation for more details.
#' @param trace whether to show the progress bar. The user is expected to have
#' set up appropriate handlers for this using the \dQuote{progressr} package.
#' @param solver choice of solver to use for the estimation of the paths.
#' @param autodiff whether to use automatic differentiation for estimation.
#' This makes use of the tsvetsad package.
#' @param ... not currently used.
#' @note The function can use parallel functionality as long as the user has set
#' up a \code{\link[future]{plan}} using the future package.
#' @return An object of class \dQuote{tsvets.profile}.
#' @aliases tsprofile
#' @method tsprofile tsvets.estimate
#' @rdname tsprofile
#' @export
#'
tsprofile.tsvets.estimate <- function(object, h = 1, nsim = 100, seed = NULL, trace = FALSE, solver = "nlminb", autodiff = FALSE, ...)
{
    sim <- simulate(object, seed = seed, nsim = nsim, h = NROW(object$spec$target$y_orig) + h)
    profile <- profile_fun(sim, object, h, trace = trace, solver = solver, autodiff = autodiff)
    return(profile)
}

profile_fun <- function(sim, object, h, cores, trace, solver, autodiff)
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
        mod <- try(estimate(spec, solver = solver, autodiff = autodiff), silent = TRUE)
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
