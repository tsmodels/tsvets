ss_vets_init <- function(y, model, damped = FALSE, frequency = 1)
{
    init_states %<-% future_lapply(1:ncol(y), function(i) {
        mod <- estimate(ets_modelspec(y[,i], model, damped = damped, normalized_seasonality = TRUE, seasonal_init = "estimate", 
                                      frequency = frequency))
        states <- matrix(mod$model$states[1,], nrow = 1)
        return(states)
    }, future.packages = c("xts","tsets"))
    init_states <- eval(init_states)
    init_states <- do.call(rbind, init_states)
    return(t(init_states))
}