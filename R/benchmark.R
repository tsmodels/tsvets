#' Model Benchmarking
#'
#' @description Benchmarks and saves the details of a model for benchmarking of
#' speed and unit testing.
#' @param object an object of class \dQuote{tsvets.spec}.
#' @param solver solver to use for estimation.
#' @param control additional controls passed to the solver.
#' @param autodiff whether to use automatic differentiation for estimation.
#' This makes use of the tsvetsad package.
#' @param ... not currently used.
#' @aliases tsbenchmark
#' @method tsbenchmark tsvets.spec
#' @rdname tsbenchmark
#' @export
#'
#'
tsbenchmark.tsvets.spec <- function(object, solver = "nlminb", control = list(trace = 0), autodiff = FALSE, ...)
{
    start <- Sys.time()
    mod <- estimate(object, solver = solver, control = control)
    end <- Sys.time()
    return(data.table(start = start, end = end, spec = list(object), estimate = list(mod), solver = solver, control = control, loglik = as.numeric(logLik(mod))))
}