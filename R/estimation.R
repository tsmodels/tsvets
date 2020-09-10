estimate.tsvets.spec = function(object, solver = "nlminb", control = list(trace = 1, iter.max = 200, eval.max = 1000), ...)
{
  env <- object$vets_env
  env$loglik <- 1
  pars <- env$pars
  pars <- pars + 1e-2
  lb <- env$lower_index
  ub <- env$upper_index
  # diagonal 2 is experimental for investigating the speed up from using a decomposed form
  # of the likelihood
  likfun <- switch(object$dependence$type,
                   "equicorrelation" = vets_llh_equicor,
                   "diagonal" = vets_llh_diag,
                   "full" = vets_llh_full,
                   "shrinkage" = vets_llh_shrink)
  tic <- Sys.time()
  # Need:
  # 1. Scaling for regressors
  # 2. Parameter transformation for unconstrained optimization
  if (solver == "optim") {
    sol <- optim(par = pars, fn = likfun, method = "L-BFGS-B", lower = lb, upper = ub, control = control, env = env)
    pars <- sol$par
    llh <- sol$value
  } else if (solver == "nlminb") {
    sol <- nlminb(start = pars, objective = likfun, control = control, lower = lb, upper = ub, env = env)
    pars <- sol$par
    llh <- sol$objective
  } else if (solver == "solnp") {
    sol <- solnp(pars = pars, fun = likfun, control = control, LB = lb, UB = ub, env = env)
    pars <- sol$pars
    llh <- tail(sol$values, 1)
  } else {
    stop("\nunrecognized solver")
  }
  sol$elapsed <- difftime(Sys.time(), tic, "minutes")
  sol$negative_llh <- llh
  filt <- vets_filter(pars, env)
  if (any(object$transform$lambda != 1)) {
    fit <- object$transform$inverse(filt$fitted[-1,], object$transform$lambda)
    fit <- xts(fit, object$target$index)
  } else {
    fit <- xts(filt$fitted[-1,], object$target$index)
  }
  object$vets_env$Amat <- filt$Amat
  object$vets_env$Fmat <- filt$Fmat
  object$vets_env$Gmat <- filt$Gmat
  out <- list(fitted = fit, States = filt$States, Error = filt$Error[-1,], spec = object, opt = sol)
  class(out) <- "tsvets.estimate"
  return(out)
}

vets_llh_equicor <- function(pars, env)
{
  # pars <- 1:8
  env$Amat[env$matrix_index] <- pars[env$parameter_index]
  Y <- env$ymat
  Amat <- t(env$Amat[,env$Amat_index[1]:env$Amat_index[2], drop = FALSE])
  Fmat <- env$Fmat
  States <- env$States
  if (env$Phi_index[1] > 0) {
    Fmat[which(as.logical(is.na(Fmat)))] <- env$Amat[,env$Phi_index[1]:env$Phi_index[2]]
  }
  if (env$X_index[1] >= 0) {
    beta <- env$Amat[,env$X_index[1]:env$X_index[2], drop = FALSE]
    xreg <- env$xreg
  } else {
    beta <- matrix(1, ncol = env$model[2], nrow = 1)
    xreg <- env$xreg
  }
  Gmat <- env$Gmat
  Hmat <- env$Hmat
  model <- env$model
  model[length(model)] <- pars[length(pars)]
  f <- vets_cpp_llh_equicor(model, as.matrix(Amat), Fmat, Hmat, Gmat, States, Y, xreg, beta)
  if (!is.finite(f$loglik) | f$condition == 1 ) {
    llh <- get("loglik", env) + 0.1*(abs(get("loglik", env)))
  } else{
    llh <- f$loglik
  }
  assign("loglik", llh, envir = env)
  return(llh)
}

vets_llh_shrink <- function(pars, env)
{
  # pars <- 1:8
  env$Amat[env$matrix_index] <- pars[env$parameter_index]
  Y <- env$ymat
  Amat <- t(env$Amat[,env$Amat_index[1]:env$Amat_index[2], drop = FALSE])
  Fmat <- env$Fmat
  States <- env$States
  if (env$Phi_index[1] > 0) {
    Fmat[which(as.logical(is.na(Fmat)))] <- env$Amat[,env$Phi_index[1]:env$Phi_index[2]]
  }
  if (env$X_index[1] >= 0) {
    beta <- env$Amat[,env$X_index[1]:env$X_index[2], drop = FALSE]
    xreg <- env$xreg
  } else {
    beta <- matrix(1, ncol = env$model[2], nrow = 1)
    xreg <- env$xreg
  }
  Gmat <- env$Gmat
  Hmat <- env$Hmat
  model <- env$model
  model[length(model)] <- pars[length(pars)]
  f <- vets_cpp_llh_shrink(model, as.matrix(Amat), Fmat, Hmat, Gmat, States, Y, xreg, beta)
  if (!is.finite(f$loglik) | f$condition == 1 ) {
    llh <- get("loglik", env) + 0.1*(abs(get("loglik", env)))
  } else{
    llh <- f$loglik
  }
  assign("loglik", llh, envir = env)
  return(llh)
}


vets_llh_diag <- function(pars, env)
{
  # pars <- 1:8
  env$Amat[env$matrix_index] <- pars[env$parameter_index]
  Y <- env$ymat
  Amat <- t(env$Amat[,env$Amat_index[1]:env$Amat_index[2], drop = FALSE])
  Fmat <- env$Fmat
  States <- env$States
  if (env$Phi_index[1] > 0) {
    Fmat[which(as.logical(is.na(Fmat)))] <- env$Amat[,env$Phi_index[1]:env$Phi_index[2]]
  }
  if (env$X_index[1] >= 0) {
    beta <- env$Amat[,env$X_index[1]:env$X_index[2], drop = FALSE]
    xreg <- env$xreg
  } else {
    beta <- matrix(1, ncol = env$model[2], nrow = 1)
    xreg <- env$xreg
  }
  Gmat <- env$Gmat
  Hmat <- env$Hmat
  model <- env$model
  f <- vets_cpp_llh_diagonal(model, as.matrix(Amat), Fmat, Hmat, Gmat, States, Y, xreg, beta)
  if (!is.finite(f$loglik) | f$condition == 1 ) {
    llh <- get("loglik", env) + 0.2*(abs(get("loglik", env)))
  } else{
    llh <- f$loglik
  }
  assign("loglik", llh, envir = env)
  return(llh)
}

vets_llh_full <- function(pars, env)
{
  # pars <- 1:8
  env$Amat[env$matrix_index] <- pars[env$parameter_index]
  Y <- env$ymat
  Amat <- t(env$Amat[,env$Amat_index[1]:env$Amat_index[2], drop = FALSE])
  Fmat <- env$Fmat
  States <- env$States
  if (env$Phi_index[1] > 0) {
    Fmat[which(as.logical(is.na(Fmat)))] <- env$Amat[,env$Phi_index[1]:env$Phi_index[2]]
  }
  if (env$X_index[1] >= 0) {
    beta <- env$Amat[,env$X_index[1]:env$X_index[2], drop = FALSE]
    xreg <- env$xreg
  } else {
    beta <- matrix(1, ncol = env$model[2], nrow = 1)
    xreg <- env$xreg
  }
  Gmat <- env$Gmat
  Hmat <- env$Hmat
  model <- env$model
  f <- vets_cpp_llh_full(model, as.matrix(Amat), Fmat, Hmat, Gmat, States, Y, xreg, beta)
  if (!is.finite(f$loglik) | f$condition == 1 ) {
    llh <- get("loglik", env) + 0.2*(abs(get("loglik", env)))
  } else{
    llh <- f$loglik
  }
  assign("loglik", llh, envir = env)
  return(llh)
}

vets_filter <- function(pars, env)
{
  env$Amat[env$matrix_index] <- pars[env$parameter_index]
  Y <- env$ymat
  Amat <- t(env$Amat[,env$Amat_index[1]:env$Amat_index[2], drop = FALSE])
  Fmat <- env$Fmat
  States <- env$States
  if (env$Phi_index[1] > 0) {
    Fmat[which(as.logical(is.na(Fmat)))] <- env$Amat[,env$Phi_index[1]:env$Phi_index[2]]
  }
  if (env$X_index[1] >= 0) {
    beta <- env$Amat[,env$X_index[1]:env$X_index[2], drop = FALSE]
    xreg <- env$xreg
  } else {
    beta <- matrix(1, ncol = env$model[2], nrow = 1)
    xreg <- env$xreg
  }
  Gmat <- env$Gmat
  Hmat <- env$Hmat
  model <- env$model
  f <- vets_cpp_filter(model, as.matrix(Amat), Fmat, Hmat, Gmat, States, Y, xreg, beta)
  f$Amat <- env$Amat
  f$Fmat <- Fmat
  f$Gmat <- Gmat
  return(f)
}
