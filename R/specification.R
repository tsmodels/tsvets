#' Model Specification
#'
#' @description Specifies an vector ETS model prior to estimation.
#' @details The specification allows to specify a vector additive damped ETS 
#' model with options for the dynamics of the states and dependence.
#' @param y an xts matrix.
#' @param level dynamics for the level component.
#' @param slope dynamics for the slope component.
#' @param damped dynamics for the dampening component.
#' @param seasonal dynamics for the seasonal component.
#' @param group a vector of indices denoting which group the series belongs to 
#' (when using the grouped dynamics).
#' @param xreg an xts matrix of external regressors.
#' @param xreg_include a matrix of dimension ncol(y) by ncol(xreg) populated with 
#' either 0, 1 or 2+ (0 = no beta, 1 = individual beta and 2 = grouped beta). 
#' It is also possible to have group wise pooling. For instance 2 variables 
#' sharing one pooled estimates, and 3 other variables sharing another grouped 
#' estimate would have values of (2,2,3,3,3). The index for group wise pooling 
#' starts at 2 and should be incremented for each new group added.
#' @param transformation a valid transformation for y from the \dQuote{tstransform} 
#' function in the \dQuote{tsaux} package (currently box-cox or logit are available).
#' @param lambda the Box Cox power transformation vector (see \code{box_cox}) 
#' in the \bold{tsaux} package. If a single NA, then it will calculate optimal 
#' lambda based on the multivariate Box Cox approach of Velila (1993), else if a 
#' vector of NA values it will calculate the individual Box Cox optimal parameters. 
#' Can also be either a single value (common lambda) or vector of values 
#' (individual lambda).
#' @param lower lower bound for the transformation.
#' @param upper upper bound for the transformation.
#' @param frequency seasonal frequency of the series.
#' @param dependence dependence structure to impose.
#' @return An object of class \dQuote{tsvets.spec} with the following slots:\cr
#' \item{target}{A list with original data series, the data series index and 
#' the sampling frequency}
#' \item{transform}{A list with details on the transformation}
#' \item{model}{A list with details the type of model dynamics}
#' \item{dependence}{A list with details about the dependence structure}
#' \item{xreg}{A list with details on the external regressors}
#' \item{vets_env}{An environment with pre-calculated state matrices and other 
#' parameters which will be passed to the estimation routine}
#' @references Athanasopoulos, G and de Silva, A. (2012),\emph{Multivariate 
#' Exponential Smoothing for Forecasting Tourist Arrivals}, Journal of Travel 
#' Research 51(5) 640–-652.\cr
#' de Silva, A., R. Hyndman, and R. D. Snyder. (2010).\emph{The Vector Innovations 
#' Structural Time Series Framework: A Simple Approach to Multivariate Forecasting}, 
#' Statistical Modelling (10) 353--74.
#' @aliases vets_modelspec
#' @rdname vets_modelspec
#' @export
#'
#'
vets_modelspec = function(y, level = c("constant","diagonal","common","full","grouped"), 
                          slope = c("none","constant","common","diagonal","full","grouped"),
                    damped = c("none","common","diagonal","full","grouped"), 
                    seasonal = c("none","common","diagonal","full","grouped"),
                    group = NULL, xreg = NULL, xreg_include = NULL, frequency = 1, 
                    transformation = "box-cox", lambda = NULL, lower = 0, upper = 1, 
                    dependence = c("diagonal","full","equicorrelation","shrinkage"))
{
  # for xreg_include  this should be a matrix with rows (i) representing y and columns (j) the regressors and populated with
  # 1s or 0s depending on whether to apply a regressor (xreg_j) to y_i. Additionally, these can be numbered so that common
  # coefficients can be included.
  if (!is.xts(y)) {
    stop("y must be an xts object")
  }
  n <- NCOL(y)
  #if (n == 1) stop("\ncannot specify a vector model with only one series")
  ynames <- colnames(y)
  if (is.null(ynames)) {
    colnames(y) <- paste0("V", 1:n)
    ynames <- colnames(y)
  }
  if (is.null(frequency)) {
    if (seasonal != "none") stop("\nseasonal models require frequency to be set.")
    frequency <- 1
  } else {
    frequency <- frequency[1]
  }
  y_orig <- y
  period <- sampling_frequency(index(y))
  if (length(transformation) != 1) {
    if (!all(transformation == "box-cox") & !all(transformation == "logit")) {
      stop("\ntransformation must be the same for all variables")
    }
  }
  transformation <- match.arg(transformation[1], c("logit","box-cox"))
  if (transformation == "logit") {
    lambda <- rep(0, n)
  } else {
    if (!is.null(lambda)) {
      if (length(lambda) != 1 & length(lambda) != n) stop("\nlambda must be of length 1 or equal to the ncol of y")
      if (all(!is.na(lambda)) & all(lambda == 1)) {
        lambda <- NULL
      }
    }
  }
  transformation_model <- list(model = transformation, lower = lower, upper = upper)
  if (!is.null(lambda)) {
    if (transformation == "box-cox") {
      if (length(lambda) == 1 && is.na(lambda)) {
        tmp <- tstransform(method = transformation, lambda = lambda, lower = lower, upper = upper, frequency = frequency, multivariate = TRUE)
        y <- tmp$transform(y)
        lambdas <- unname(attr(y, "lambda"))
        transform <-  lapply(1:n, function(i) tstransform(method = transformation, lambda = lambdas[i], frequency = frequency, lower = lower, upper = upper))
        for (i in 1:n) {
          transform[[i]]$lambda <- lambdas[i]
          transform[[i]]$name <- "box-cox"
        }
      } else if (length(lambda) == n) {
        transform <-  lapply(1:n, function(i) tstransform(method = transformation, lambda = lambda[i], frequency = frequency, lower = lower, upper = upper))
        for (i in 1:n) {
            tmp <- transform[[i]]$transform(y[,i], lambda = lambda[i], frequency = frequency, lower = lower, upper = upper)
            transform[[i]]$lambda <- unname(attr(tmp,"lambda"))
            transform[[i]]$name <- "box-cox"
            y[,i] <- as.numeric(tmp)
        }
      } else if (length(lambda) == 1 & !is.na(lambda)) {
        transform <-  lapply(1:n, function(i) tstransform(method = transformation, lambda = lambda, frequency = frequency, lower = lower, upper = upper))
        for (i in 1:n) {
          tmp <- transform[[i]]$transform(y[,i], lambda = lambda, frequency = frequency, lower = lower, upper = upper)
          transform[[i]]$lambda <- unname(attr(tmp,"lambda"))
          transform[[i]]$name <- "box-cox"
          y[,i] <- as.numeric(tmp)
        }
      }
    } else {
      transform <-  lapply(1:n, function(i) tstransform(method = transformation, lower = lower, upper = upper))
      for (i in 1:n) {
        transform[[i]]$lambda <- 0
        tmp <- transform[[i]]$transform(y[,i])
        transform[[i]]$name <- "logit"
        y[,i] <- as.numeric(tmp)
      }
    }
  } else {
    transform <- NULL
  }
  # check and index missing values
  good_matrix <- matrix(1, ncol = ncol(y), nrow = nrow(y))
  good_index <- rep(0, nrow(y))
  if (any(is.na(y))) {
    exc <- which(is.na(y), arr.ind = TRUE)
    good_matrix[exc] <- NA
    nm <- NROW(na.omit(good_matrix))
    if (nm <= ncol(good_matrix)) stop("\nTotal common non-missingness is less than number of variables. Not currently supported.")
    good_matrix <- na.fill(good_matrix, fill = 0)
  }
  good_index[which(rowSums(good_matrix) == ncol(y))] <- 1
  
  level <- match.arg(arg = level[1], choices = c("constant","diagonal","common","full","grouped"))
  slope <- match.arg(arg = slope[1], choices = c("none","constant","common","diagonal","full","grouped"))
  damped <- match.arg(arg = damped[1], choices = c("none","common","diagonal","full","grouped"))
  seasonal <- match.arg(arg = seasonal[1], choices = c("none","common","diagonal","full","grouped"))
  dependence <- match.arg(arg = dependence[1], choices = c("diagonal","full","equicorrelation","grouped_equicorrelation","shrinkage"))
  if (dependence == "grouped_equicorrelation") stop("\ngrouped_equicorrelation not yet implemented")
  if (any(c(level,slope,damped,seasonal) %in% "grouped")) {
    if (is.null(group)) stop("\ngroup cannot be NULL for grouped choice.")
    if (length(group) != n) stop("\nlength of group vector must be equal to number of cols of y.")
    if (max(group) > n) stop("\nmax group > ncol y...check and resubmit.")
    if (all(group == max(group))) stop("\ngroup variable is the same for all y. Try common instead.")
  } else {
    group <- NULL
  }
  L <- list(frequency = frequency, n = ncol(y), group = group, slope = ifelse(slope != "none", TRUE, FALSE), 
            seasonal = ifelse(seasonal != "none", TRUE, FALSE), series_names = ynames, xreg = xreg, xreg_include = xreg_include)
  # H matrix
  Hmat <- ss_mat_H(L)
  n_states <- nrow(Hmat)
  # G matrix
  Gmat <- ss_mat_G(L)
  # A matrix
  k <- 0
  model <- "ANN"
  A <- ss_mat_A(L, level, counter = k)
  AA <- rbind(Matrix(A$parameters, L$n, L$n))
  i_index <- A$index
  e_index <- A$estimate
  lower_index <- A$lower
  upper_index <- A$upper
  parnames <- A$parnames
  if (!all(is.na(A$index))) {
    k <- max(A$index, na.rm = TRUE)
  }
  # B matrix
  if (slope != "none") {
    substr(model,2,2) <- "A"
    B <- ss_mat_B(L, slope, counter = k)
    if (!all(is.na(B$index))) {
      k <- max(B$index, na.rm = TRUE)
    }
    AA <- rbind(AA, Matrix(B$parameters, L$n, L$n))
    i_index <- c(i_index, B$index)
    e_index <- c(e_index, B$estimate)
    lower_index <- c(lower_index, B$lower)
    upper_index <- c(upper_index, B$upper)
    parnames <- c(parnames, B$parnames)
    
  } else {
    B <- NULL
  }
  # K matrix
  if (seasonal != "none") {
    substr(model,3,3) <- "A"
    K <- ss_mat_K(L, seasonal, counter = k)
    AA <- rbind(AA, Matrix(K$parameters, L$frequency * L$n, L$n, byrow = TRUE))
    # AA <- t(AA)
    # AA[which(is.na(AA))] <- pars[i_index[which(i_index!=0)]]
    i_index <- c(i_index, K$index, rep(K$index, L$frequency - 1))
    e_index <- c(e_index, K$estimate, rep(0, L$n * (L$frequency - 1)))
    lower_index <- c(lower_index, K$lower, rep(0, L$n * (L$frequency - 1)))
    upper_index <- c(upper_index, K$upper, rep(0, L$n * (L$frequency - 1)))
    parnames <- c(parnames, K$parnames, rep(as.numeric(NA), L$n * (L$frequency - 1)))
    if (!all(is.na(B$index))) {
      k <- max(K$index, na.rm = TRUE)
    }
  }
  Amat_index <- c(1, nrow(AA))
  # Phi Matrix
  if (slope != "none") {
    Phi <- ss_mat_Phi(L, damped, counter = k)
    
    if (damped != "none") {
      Phi_index <- c(nrow(AA) + 1, nrow(AA) + L$n)
    } else {
      Phi_index <- c(-1, -1)
    }
    AA <- rbind(AA, Matrix(Phi$parameters, L$n, L$n, byrow = TRUE))
    i_index <- c(i_index, Phi$index)
    e_index <- c(e_index, Phi$estimate)
    lower_index <- c(lower_index, Phi$lower)
    upper_index <- c(upper_index, Phi$upper)
    parnames <- c(parnames, Phi$parnames)
    k <- max(i_index, na.rm = TRUE)
  } else {
    Phi <- NULL
    Phi_index <- c(-1, -1)
  }
  
  # Initial states [currently we do not allow estimation of these as it explodes the number of parameters]
  init_states <- ss_vets_init(y, model, damped = ifelse(damped == "none", FALSE, TRUE), frequency = frequency)
  # Remove indices now from y since they are no longer needed
  y <- coredata(y)
  y_index <- index(y_orig)
  y_orig <- coredata(y_orig)
  
  S_index <- c(nrow(AA) + 1, nrow(AA) + 1)
  AA <- rbind(AA, Matrix(init_states[1,,drop = FALSE]))
  
  i_index <- c(i_index, rep(0, n))
  e_index <- c(e_index, rep(0, n))
  lower_index <- c(lower_index, rep(0, n))
  upper_index <- c(upper_index, rep(0, n))
  parnames <- c(parnames, paste0("l0[",ynames,"]"))
  
  if (slope != "none") {
      S_index[2] <- S_index[2] + 1
      AA <- rbind(AA, Matrix(init_states[2,,drop = FALSE]))
      i_index <- c(i_index, rep(0, n))
      e_index <- c(e_index, rep(0, n))
      lower_index <- c(lower_index, rep(0, n))
      upper_index <- c(upper_index, rep(0, n))
      parnames <- c(parnames, paste0("b0[",ynames,"]"))
  }
  
  if (seasonal != "none") {
      S_index[2] <- S_index[2] + frequency
      AA <- rbind(AA, Matrix(tail(init_states, frequency)))
      i_index <- c(i_index, rep(0, n * frequency))
      e_index <- c(e_index, rep(0, n * frequency))
      lower_index <- c(lower_index, rep(0, n * frequency))
      upper_index <- c(upper_index, rep(0, n * frequency))
      parnames <- c(parnames, paste0("s",rep(0:(frequency - 1), each = n),"[",ynames,"]"))
  }
  # xreg
  if (!is.null(xreg)) {
    # Check indices of xreg against y
    xreg <- check_xreg(xreg, y_index)
    xreg <- coredata(xreg)
    if (is.null(colnames(xreg))) colnames(xreg) <- paste0("X", 1:ncol(xreg))
    L$xreg_include <- xreg_include
    L$xreg <- xreg
    X <- ss_mat_xreg(L, counter = k)
    beta <- Matrix(as.numeric(X$parameters), ncol = n, nrow = ncol(xreg), byrow = TRUE)
    i_index <- c(i_index, X$index)
    e_index <- c(e_index, X$estimate)
    lower_index <- c(lower_index, X$lower)
    upper_index <- c(upper_index, X$upper)
    parnames <- c(parnames, paste0("X[", rep(colnames(xreg), each = n),"][",ynames,"]"))
    X_index <- c(nrow(AA) + 1, nrow(AA) + NROW(beta))
    AA <- rbind(AA,  beta)
    include_xreg <- TRUE
  } else {
    X_index <- c(-1, -1)
    xreg <- matrix(0, ncol = 1, nrow = NROW(y))
    beta <- matrix(0, nrow = ncol(y), ncol = 1)
    include_xreg <- FALSE
  }
  # Fmat
  if (damped != "none") { 
    Phimat <- matrix(NA, L$n, L$n)
  } else {
    Phimat <- diag(1, L$n, L$n)
  }
  Fmat <- ss_mat_F(L, Phimat)
  if (damped != "none") {
    Fmat_index <- which(is.na(as.vector(Fmat)))
  } else {
    Fmat_index <- c(-1, -1)
  }
  # dependence (we split out the functions depending on the dependence type as it is simpler to
  # enhance and maintain) without affecting the other code
  if (dependence == "equicorrelation") {
    dL <- list(R = 0.5, lower = -0.99, upper = 0.99, type = dependence, estimate = TRUE)
  } else if (dependence == "shrinkage") {
    dL <- list(R = 0, lower = 0, upper = 1, type = dependence, estimate = TRUE)
  } else{
    dL <- list(R = as.numeric(NA), lower = as.numeric(NA), upper = as.numeric(NA), type = dependence, estimate = FALSE)
  }
  rownames(AA) <- NULL
  vets_env <- new.env()
  vets_env$ymat <- t(rbind(matrix(0, nrow = 1, ncol = ncol(y)), coredata(y)))
  vets_env$Amat <- as.matrix(t(AA))
  vets_env$Hmat <- t(Hmat)
  vets_env$Fmat <- Fmat
  vets_env$Gmat <- Gmat
  vets_env$Amat_index <- Amat_index
  vets_env$Fmat_index <- Fmat_index
  vets_env$Phi_index <- Phi_index
  vets_env$S_index <- S_index
  vets_env$X_index <- X_index
  vets_env$xreg <- t(rbind(matrix(0, nrow = 1, ncol = ncol(xreg)), coredata(xreg)))
  vets_env$parameter_index <- i_index[which(i_index > 0)]
  vets_env$matrix_index <- which(is.na(as.vector(t(AA))))
  vets_env$estimate_index <- e_index
  vets_env$lower_index <- lower_index[which(e_index == 1)]
  vets_env$upper_index <- upper_index[which(e_index == 1)]
  vets_env$pars <- rep(0, length(which(e_index == 1)))
  vets_env$parnames <- parnames[which(e_index ==  1)]
  if (dependence == "equicorrelation") {
    vets_env$pars <- c(vets_env$pars, 0.5)
    vets_env$parnames <- c(vets_env$parnames, "rho")
    vets_env$lower_index <- c(vets_env$lower_index, 0)
    vets_env$upper_index <- c(vets_env$upper_index,  0.99)
  } else if (dependence == "shrinkage") {
    vets_env$pars <- c(vets_env$pars, 0)
    vets_env$parnames <- c(vets_env$parnames, "rho")
    vets_env$lower_index <- c(vets_env$lower_index, 1e-12)
    vets_env$upper_index <- c(vets_env$upper_index, 1)
  }
  States <- matrix(as.vector(vets_env$Amat[,vets_env$S_index[1]:vets_env$S_index[2]]), ncol = 1)
  States <- cbind(States, matrix(0, ncol = NROW(y), nrow = nrow(States)))
  vets_env$States <- States
  vets_env$model <- c(nrow(y) + 1, ncol(y), as.integer(include_xreg), 0.5)
  vets_env$init_states <- init_states
  # vets.env$Amat[vets.env$matrix.index] <- p[vets.env$parameter.index]
  # extract
  spec <- list()
  spec$target$y <- y
  spec$target$y_orig <- y_orig
  spec$target$index <- y_index
  spec$target$frequency <- frequency
  spec$target$sampling <- period
  spec$target$y_names <- ynames
  spec$target$good_matrix <- good_matrix
  spec$target$good_index <- good_index
  spec$transform <- transform
  spec$model$level <- level
  spec$model$slope <- slope
  spec$model$damped <- damped
  spec$model$seasonal <- seasonal
  spec$model$group <- group
  spec$dependence <- dL
  spec$xreg$xreg <- xreg
  spec$xreg$xreg_include <- xreg_include
  spec$xreg$include_xreg <- include_xreg
  spec$vets_env <- vets_env
  class(spec) <- "tsvets.spec"
  return(spec)
}


#' Model Specification Extractor
#'
#' @description Extracts a model specification (class \dQuote{tsvets.spec}) from
#' an object of class \dQuote{tsvets.estimate}.
#' @param object an object of class \dQuote{tsvets.estimate}.
#' @param y an optional new xts vector.
#' @param lambda an optional lambda parameter for the Box Cox transformation (if
#' previously used).
#' @param xreg an optional matrix of regressors.
#' @param ... not currently used.
#' @note This function is used by other functions in the package such as the
#' backtest which requires rebuilding the specification for each re-estimation
#' step with updated data but keeping all else equal.
#' @return An object of class \dQuote{tsvets.spec}.
#' @aliases tsspec
#' @method tsspec tsvets.estimate
#' @rdname tsspec
#' @export
#'
#'
tsspec.tsvets.estimate <- function(object, y = NULL, lambda = NULL, xreg = NULL, ...)
{
  if (is.null(y)) {
    y <- xts(object$spec$target$y_orig, object$spec$target$index)
  }
  if (!is.null(xreg)) {
    xreg <- coredata(xreg)
    if (nrow(xreg) != NROW(y)) stop("\nxreg should have the same number of rows as y")
    if (ncol(xreg) > (NROW(y)/2)) warning("\nnumber of regressors greater than half the length of y!")
  } else {
    if (object$spec$xreg$include_xreg) {
      xreg <- object$spec$xreg$xreg
      if (nrow(xreg) != NROW(y)) stop("\nxreg should have the same number of rows as y")
      if (ncol(xreg) > (NROW(y)/2)) warning("\nnumber of regressors greater than half the length of y!")
    } else {
      xreg <- NULL
    }
  }
  if (is.null(lambda)) {
    lambda <- sapply(1:length(object$spec$transform), function(i) object$spec$transform[[i]]$lambda)
  }
  spec <- vets_modelspec(y, level = object$spec$model$level, 
                 slope = object$spec$model$slope, 
                 damped = object$spec$model$damped, 
                 seasonal = object$spec$model$seasonal,
                 group = object$spec$model$group, 
                 xreg = xreg, 
                 xreg_include = object$spec$xreg$xreg_include, 
                 frequency = object$spec$target$frequency, 
                 transformation = object$spec$transform[[1]]$name,
                 lambda = lambda, lower = object$spec$transform[[1]]$lower, 
                 upper = object$spec$transform[[1]]$upper, 
                 dependence = object$spec$dependence$type)
  return(spec)
}

