estimate.tsvets.spec = function(object, solver = "nlminb", control = list(trace = 0, iter.max = 200, eval.max = 1000), autodiff = FALSE, ...)
{
  env <- object$vets_env
  env$loglik <- 1
  pars <- env$pars
  pars <- pars + 1e-2
  lb <- env$lower_index
  ub <- env$upper_index
  env$good <- t(object$target$good_matrix)
  env$select <- t(object$target$good_index)
  if (any(is.na(env$ymat))) {
    env$ymat <- na.fill(env$ymat, fill = 0)
  }
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
  hessian <- NULL
  gradient <- NULL
  if (autodiff) {
    solver <- match.arg(solver[1], choices = c("nlminb","optim"))
    opt <- estimate_ad(object, solver = solver, control = control, ...)
    pars <- opt$pars
    hessian <- opt$hessian
    gradient <- opt$gradient
    sol <- opt$sol
    llh <- opt$llh
  } else {
    solver <- match.arg(solver[1], choices = c("nlminb","optim","solnp","gosolnp","nloptr"))
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
    } else if (solver == "gosolnp") {
        sol <- gosolnp(pars = pars, fun = likfun, control = control, LB = lb, UB = ub, env = env, ...)
        pars <- sol$pars
        llh <- tail(sol$values, 1)
    } else if (solver == "nloptr") {
      sol <- nloptr(x0 = pars, eval_f = likfun, lb = lb, ub = ub, env = env, opts = list("algorithm" = "NLOPT_GN_MLSL_LDS", xtol_rel = 1e-8, maxeval = 1000, maxtime = 60*2,
                                                                                         local_opts = list(algorithm = "NLOPT_LN_NELDERMEAD"), print_level = 0))
      pars <- sol$solution
      if (is.null(control$maxeval)) maxeval <- 2000 else maxeval <- control$maxeval
      if (is.null(control$xtol_rel)) xtol_rel <- 1e-8 else xtol_rel <- control$xtol_rel
      sol <- nloptr(x0 = pars, eval_f = likfun, lb = lb, ub = ub, env = env, opts = list("algorithm" = "NLOPT_LN_COBYLA", maxeval = maxeval, xtol_rel = xtol_rel, print_level = as.integer(control$trace)))
      pars <- sol$solution
      llh <- sol$objective
      sol$par <- pars
    } else {
      stop("\nunrecognized solver")
    }
  }
  sol$elapsed <- difftime(Sys.time(), tic, "minutes")
  sol$negative_llh <- llh
  filt <- vets_filter(pars, env)
  if (!is.null(object$transform)) {
      # need the minverse
      fit <- do.call(cbind, lapply(1:length(object$transform), function(i) object$transform[[i]]$inverse(filt$fitted[-1,i], object$transform[[i]]$lambda)))
  } else {
      fit <- filt$fitted[-1,]
  }
  fit <- xts(fit, object$target$index)
  
  object$vets_env$Amat <- filt$Amat
  object$vets_env$Fmat <- filt$Fmat
  object$vets_env$Gmat <- filt$Gmat
  out <- list(fitted = fit, States = filt$States, Error = filt$Error[-1,], gradient = gradient, hessian = hessian, spec = object, opt = sol, autodiff = autodiff)
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
  good <- env$good
  select <- as.numeric(env$select)
  f <- vets_cpp_llh_equicor(model, as.matrix(Amat), Fmat, Hmat, Gmat, States, Y, xreg, beta, good, select)
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
  good <- env$good
  select <- as.numeric(env$select)
  model[length(model)] <- pars[length(pars)]
  f <- vets_cpp_llh_shrink(model, as.matrix(Amat), Fmat, Hmat, Gmat, States, Y, xreg, beta, good, select)
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
  good <- env$good
  select <- as.numeric(env$select)
  model <- env$model
  f <- vets_cpp_llh_diagonal(model, as.matrix(Amat), Fmat, Hmat, Gmat, States, Y, xreg, beta, good, select)
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
  good <- env$good
  select <- as.numeric(env$select)
  model <- env$model
  f <- vets_cpp_llh_full(model, as.matrix(Amat), Fmat, Hmat, Gmat, States, Y, xreg, beta, good, select)
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
  good <- env$good
  model <- env$model
  f <- vets_cpp_filter(model, as.matrix(Amat), Fmat, Hmat, Gmat, States, Y, xreg, beta, good)
  f$Amat <- env$Amat
  f$Fmat <- Fmat
  f$Gmat <- Gmat
  return(f)
}
