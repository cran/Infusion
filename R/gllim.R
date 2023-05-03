.loggausspdf <- function (X, mu, solve_t_chol_sigma) {
  if (ncol(mu) == ncol(X)) {
    X = X - mu
  } else {
    X = sweep(X, 1, mu, "-")
  }
  Q <- solve_t_chol_sigma %*% X
  q <-  colSums(Q^2)
  # c = d * log(2 * pi) + 2 * sum(log(diag(U)))
  c <- nrow(X) * log(2 * pi) - 2 *attr(solve_t_chol_sigma,"logdet")
  -(c + q)/2
}

.logsumexp <- function (x) {
  y = apply(x, 1, max)
  x = x - y
  s = y + log(rowSums(exp(x)))
  i = is.infinite(y)
  if (sum(i) > 0) s[i] = y[i]
  return(s)
}

.wrap_gllim <- function(RGPpars, summstats, nbCluster, cstr=list(Sigma="i")) {
  if (length(nbCluster)>1L) {
    objfn <- function(RGPpars,summstats, K, cstr=list(Sigma="i")) {
      cluObject <- .do_call_wrap("gllim",arglist=list(tapp=RGPpars, yapp=summstats,in_K=K, cstr=cstr), 
                                 pack="xLLiM")
      logL <- cluObject$LLf
      df <- cluObject$nbpar
      AIC <- -2*logL+2*df
      return(list(obj=AIC, cluObject=cluObject))
    }
    resu <- .minImIze(init=NULL, objfn=objfn, range=nbCluster, trace=FALSE, RGPpars=RGPpars, 
                      summstats=summstats, cstr=cstr)$cluObject   
  } else resu <- .do_call_wrap("gllim",arglist=list(tapp=RGPpars, yapp=summstats,in_K=nbCluster, cstr=cstr), pack="xLLiM")
  
  nclu <- length(resu$pi)
  
  RGP_COVlist <- logl_COVlist <- vector("list", nclu)
  for (k in seq_len(nclu)) {
    COV <- resu$Gamma[,,k] 
    # coerced to scalar if dims are 1, while drop=FALSE would keep it as array. Neither is OK , we need matrix =>
    if (is.null(dim(COV)))  dim(COV) <- c(1L,1L)
    colnames(COV) <- rownames(COV) <- rownames(RGPpars) # fittedPars
    RGP_COVlist[[k]] <- COV
    #
    COV <- resu$Sigma[,,k] 
    if (is.null(dim(COV)))  dim(COV) <- c(1L,1L)
    logl_COVlist[[k]] <- COV
  }
  resu$RGP_COVlist <- RGP_COVlist
  resu$logl_COVlist <- logl_COVlist
  
  resu$clu_params$solve_t_chol_sigma_lists <- list(
    logldens=apply(resu$Sigma, 3, .solve_t_cholfn, simplify=FALSE), # of the N(y: Ak.x+bk, Sigmak) => log logLik of parms 'x' when used as function of x
    RGPdens=lapply(RGP_COVlist, .solve_t_cholfn) # of the N(x: ck, Gammak) => marginal distrib of params in the reftable
  )
  
  class(resu) <- c("gllim",class(resu))
  resu
}

# Algebra as in xLLiM::gllim_inverse_map which also compute the pars of the reciprocal conditional, never needed in Infusion

.logL.gllim <- function(object, # gllim result
                        newRGPpars, # Reftable Generating Process parameters. Must have dim
                        projstats, log) {
  np <- nrow(object$c) # L # number of (projected) stats
  nclu <- ncol(object$c) # K
  size <- nrow(newRGPpars) # N
  proj <- matrix(NA, nrow=nclu, ncol = size)
  condcluprobs <- matrix(0, nrow=size, ncol = nclu) # one conditional vector for each row of newRGPpars ((eta_k) in KuglerFD22)
  solve_t_chol_RGP <- object$clu_params$solve_t_chol_sigma_lists$RGPdens
  solve_t_chol_logl <- object$clu_params$solve_t_chol_sigma_lists$logldens
  for (k in 1:nclu) {
    if (np == 1L) {
      Ak = object$A[, , k, drop = FALSE]
    } else Ak = matrix(object$A[, , k], nrow = np, ncol = length(projstats))
    muk <- tcrossprod(Ak, newRGPpars) + object$b[, k]
    proj[k, ] <- .fast_dmvnorm(x = projstats, mean = muk, 
                               solve_t_chol_sigma=solve_t_chol_logl[[k]], log = FALSE)
    condcluprobs[, k] <- log(object$pi[k]) + .loggausspdf(t(newRGPpars), mu=object$c[, k, drop = FALSE], solve_t_chol_RGP[[k]])
  }
  den <- .logsumexp(condcluprobs)
  condcluprobs <- sweep(condcluprobs, 1, den, "-")
  condcluprobs <- exp(condcluprobs)
  x_exp <- proj * t(condcluprobs)
  x_exp <- colSums(x_exp)
  if (log) {
    return(log(x_exp))
  } else return(x_exp)
}

# fitted 'reftable' density of the predictors, which are the DGP parameters in the Infusion context 
# mixture of (N(., c_k, Gamma_k))
.RGPdens.gllim <- function (object, # gllim object
                            newRGPpars, log) {
  nsim <- nrow(newRGPpars) # N
  nclu <- ncol(object$c) # K
  proj <- matrix(NA, ncol = nsim, nrow=nclu)
  solve_t_chol_sigma_list <- object$clu_params$solve_t_chol_sigma_lists$RGPdens
  for (k in 1:nclu) {
    proj[k, ] <- .fast_dmvnorm(newRGPpars, mean=object$c[, k], 
                               solve_t_chol_sigma=solve_t_chol_sigma_list[[k]], log = FALSE)
  }
  x_exp <- proj * object$pi
  x_exp = colSums(x_exp)
  if (log) {
    return(log(x_exp))
  } else return(x_exp)
}

.get_dens_from_GMM.gllim <- function(X, # parameters only 
                                    object, # SLik_j object
                                    tstat.obs, # 1-row matrix as otherwise more cases should be considered for cbind'ing
                                    log, 
                                    which, # "safe" version ignores, by correcting, spuriously high logL in area of low parameter density.
                                    ...)  {
  if (inherits(X,"data.frame")) {
    X <- matrix(unlist(X,recursive = FALSE, use.names = FALSE), nrow=nrow(X),
                dimnames=list(NULL,colnames(X)) ) # newdata <- as.matrix(newdata)
  } else if (is.null(dim(X))) dim(X) <- c(1L, length(X))
  solve_t_chol_sigma_lists <- object$gllimobj$lisolve_t_chol_sigma_lists
  if (which != "parvaldens") {
    condvaldens <- .logL.gllim(object$gllimobj, newRGPpars=X, projstats=tstat.obs, 
                               log=log)
  }
  if (which=="lik") return(condvaldens)
  parvaldens <- .RGPdens.gllim(object$gllimobj, newRGPpars=X, 
                                log=log)
  if (which == "parvaldens") return(parvaldens)
  if (which=="safe") {
    if (log) {
      condvaldens <- condvaldens + pmin(0,parvaldens-object$thr_dpar) # decreasing returned logL when parvaldens<object$thr_dpar
    } else {
      condvaldens <- condvaldens*pmin(1,parvaldens/object$thr_dpar)
    }
    return(condvaldens)
  }
  if (which=="jointvaldens") { 
    if (log) {
      jointvaldens <- condvaldens + parvaldens
    } else {
      jointvaldens <- condvaldens*parvaldens
    }
    return(jointvaldens)
  }
  stop("Unhandled 'which' value")
}



..gllim.condsimul.stats <-  function (object, # gllim object
                                      RGPpars, size, cbind., colTypes) {
  if (is.null(dim(RGPpars))) {
    dim(RGPpars) <- c(1L, length(RGPpars)) 
    # input 'size' required
  } else { # parameter matrix, one sample per row
    size <- 1L
  }
  npred <- ncol(RGPpars)
  nsim <- nrow(RGPpars)
  nresp = nrow(object$c)
  nclu = ncol(object$c)
  simuls = matrix(NA, nrow = size*nsim, ncol=nresp)
  condcluprobs <- matrix(0, nrow=nsim, ncol = nclu) # one conditional vector for each row of newRGPpars ((eta_k) in KuglerFD22)
  Ak_list <- bk_list <- vector("list", nclu)
  solve_t_chol_RGP <- object$clu_params$solve_t_chol_sigma_lists$RGPdens
  for (k in 1:nclu) {
    if (nresp == 1L) {
      Ak_list[[k]] <- object$A[, , k, drop = FALSE]
    } else Ak_list[[k]] <- matrix(object$A[, , k], ncol = npred, nrow = nresp)
    bk_list[[k]] <- object$b[, k]
    condcluprobs[, k] <- log(object$pi[k]) + .loggausspdf(t(RGPpars), mu =  object$c[, k, drop = FALSE], solve_t_chol_RGP[[k]])
  }
  den <- .logsumexp(condcluprobs)
  condcluprobs = sweep(condcluprobs, 1, den, "-")
  condcluprobs = exp(condcluprobs)
  seq_nclu <- seq(nclu)
  logl_COVlist <- object$logl_COVlist
  if (nsim==1L) {
    rclu <- sample(seq_nclu, size=size, prob=condcluprobs, replace=TRUE) 
    rclutable <- table(rclu)
    simuls <- vector("list", length(rclutable))
    for (char_k in names(rclutable)) {
      k <- as.numeric(char_k)
      mu <- Ak_list[[k]] %*% RGPpars[1,] + bk_list[[k]]
      simuls[[char_k]] <- rmvnorm(rclutable[char_k], mean=mu, sigma= logl_COVlist[[k]])
    }
    simuls <- do.call(rbind, simuls)
  } else {
    for (i in seq_len(nsim))  {
      rclu <- sample(seq_nclu, size=1, prob=condcluprobs[i,]) 
      mu <- Ak_list[[rclu]] %*% RGPpars[i,] + bk_list[[rclu]]
      simuls[i,] <- rmvnorm(1L, mean=mu, sigma= logl_COVlist[[k]])
    }
  }
  
  colnames(simuls) <- colTypes$statNames
  if (cbind.) {
    colnames(RGPpars) <- colTypes$fittedPars
    simuls <- cbind(RGPpars, simuls)
  }
  return(simuls)
}

# simulation of summ stats conditional on RGPpars
.gllim.condsimul.stats <- function(object, # gllim object
                                   RGPpars,
                            seed=NULL, 
                            size=1L, # number of points for each simulation 
                            drop,
                            cbind.,
                            colTypes, 
                            ...
                            ) {
  ## RNG stuff copied from simulate.lm
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if ( ! is.null(seed)) {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  #
  if (inherits(RGPpars,"list")) {
    resu <- vector("list",  length(RGPpars))
    for (it in seq_along(resu)) {
      resu[[it]] <- ..gllim.condsimul.stats(object, RGPpars[[it]], size=size, cbind.=cbind., colTypes=colTypes) 
    }
  } else resu <- ..gllim.condsimul.stats(object, RGPpars, size=size, cbind.=cbind., colTypes=colTypes)

  resu
}

..simulate.gllimX <- function(object, size, colTypes) {
  nclu <- length(object$pi)
  rclu <- sample(seq(nclu), size=size, prob=object$pi, replace=TRUE) 
  rclutable <- table(rclu)
  simuls <- vector("list", length(rclutable))
  RGP_COVlist <- object$RGP_COVlist
  for (char_k in names(rclutable)) {
    k <- as.numeric(char_k)
    simuls[[char_k]] <- rmvnorm(rclutable[char_k], mean=object$c[, k], sigma= RGP_COVlist[[k]])
  }
  simuls <- do.call(rbind, simuls)
  colnames(simuls) <- colTypes$fittedPars
  return(simuls)
}

.simulate.gllimX <- function(object, 
                            nsim=1L, # size of list of bootstrap replicates
                            seed=NULL, 
                            size=1L, # number of points for each simulation 
                            colTypes, 
                            ...
) {
  ## RNG stuff copied from simulate.lm
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if ( ! is.null(seed)) {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  #
  if (nsim>1L) {
    resu <- vector("list",  nsim)
    for (it in seq_along(resu)) {
      resu[[it]] <- ..simulate.gllimX(object, size=size, colTypes=colTypes) 
    }
    resu
  } else ..simulate.gllimX(object, size=size, colTypes=colTypes) 
}

.condProfoutParMeans.gllim <- function(gllimobj, fittedPars, 
                                 given, expansion=Infusion.getOption("expansion")) {  # expansion=1 to get the conditional distribution
  nbCluster <- length(gllimobj$pi)
  MEAN <- t(gllimobj$c)
  givenNames <- names(given)
  colnames(MEAN) <- fittedPars
  For <- setdiff(fittedPars,givenNames) 
  condmeans <- MEAN[,For,drop=FALSE] # resizes, but will be modified
  RGP_COVlist <- gllimobj$RGP_COVlist
  for (clu_it in seq_len(nbCluster)) {
    COV <- RGP_COVlist[[clu_it]]
    sig22 <-  COV[givenNames,givenNames,drop=FALSE]
    sig12 <-  COV[For,givenNames,drop=FALSE]
    mean2 <- MEAN[clu_it,givenNames]
    condmeans[clu_it,] <- MEAN[clu_it,For] + sig12 %*% solve(sig22,given-mean2)
  }
  condmeans # matrix nclu * # profiledOutPars
}



