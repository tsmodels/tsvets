#' @rawNamespace useDynLib(tsvets,.registration = TRUE)
#' @import tsmethods
#' @import methods
#' @import data.table
#' @importFrom utils head tail data
#' @importFrom stats median na.omit fitted coef quantile residuals predict optim nlminb logLik cov cor cov2cor sd simulate var
#' @importFrom graphics grid legend lines par plot points abline axis axis.Date axis.POSIXct box polygon layout mtext
#' @importFrom grDevices gray colorRampPalette n2mfrow
#' @importFrom zoo index as.zoo zoo coredata na.locf na.fill
#' @importFrom future.apply future_lapply
#' @importFrom future %<-%
#' @importFrom progressr handlers progressor
#' @importFrom nloptr nloptr
#' @importFrom MASS mvrnorm
#' @importFrom corpcor is.positive.definite make.positive.definite
#' @importFrom tsaux mape mase bias mslre wape mis wslre do.call.fast box_cox check_xreg check_newxreg future_dates sampling_frequency tstransform
#' @importFrom viridis viridis_pal
#' @importFrom MVN mvn
#' @importFrom Rsolnp solnp gosolnp
#' @importFrom corpcor is.positive.definite make.positive.definite
#' @importFrom Matrix Matrix t
#' @importFrom xts xts as.xts is.xts endpoints
#' @importFrom tsets ets_modelspec
#' @importFrom tsvetsad estimate_ad.tsvets.spec
#' @importFrom Rcpp evalCpp loadModule
"_PACKAGE"