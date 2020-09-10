ss_vets_init <- function(y, model, damped = FALSE, frequency = 1, cores = 1)
{
    if (cores == 1) {
        i <- 1
        init_states <- foreach(i = 1:ncol(y), .combine = rbind) %do% {
            mod <- estimate(ets_modelspec(y[,i], model, damped = damped, normalized_seasonality = TRUE, seasonal_init = "estimate", 
                                          frequency = frequency))
            states <- matrix(mod$model$states[1,], nrow = 1)
            return(states)
        }
    } else {
        cl <- makeCluster(cores)
        registerDoSNOW(cl)
        clusterExport(cl,"y", envir = environment())
        clusterEvalQ(cl,"library(tsets)")
        clusterEvalQ(cl,"library(xts)")
        i <- 1
        init_states <- foreach(i = 1:ncol(y), .combine = rbind) %dopar% {
            mod <- estimate(ets_modelspec(y[,i], model, damped = damped, normalized_seasonality = TRUE, seasonal_init = "estimate", 
                                          frequency = frequency))
            states <- matrix(mod$model$states[1,], nrow = 1)
            return(states)
        }
        stopCluster(cl)
        
    }
    return(t(init_states))
}