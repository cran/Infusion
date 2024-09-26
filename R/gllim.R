.loggausspdf <- function (X, mu, solve_t_chol_sigma) {
  if (ncol(mu) == ncol(X)) {
    X = X - mu
  } else {
    X = sweep(X, 1, mu, "-")
  }
  Q <- solve_t_chol_sigma %*% X
  q <-  colSums(Q^2)
  # c = d * log(2 * pi) + 2 * sum(log(diag(U)))
  c <- nrow(X) * .log2pi - 2 *attr(solve_t_chol_sigma,"logdet")
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

.wrap_minimize_over_nbClu <- function(init=NULL, objfn, range, trace=FALSE, ...) {
  fneval <- 0L
  lower <- min(range)
  upper <- max(range)
  range <- upper - lower
  if (is.null(init)) init <- floor(lower+range/2L)
  if (trace) print(c(lower=lower,upper=upper,fneval=fneval,init=init))
  objlist <- list()
  while (range>0L) {
    vcinit <- as.character(init)
    if (is.null(objlist[[vcinit]])) {
      objlist[[vcinit]] <- objfn(init, ...)
      fneval <- fneval+1L
    }
    objinit <- objlist[[vcinit]]$obj
    #
    vp <- min(upper, init+1L)
    vc <- as.character(vp)
    if (is.null(objlist[[vc]])) {
      objlist[[vc]] <- objfn(vp, ...)
      fneval <- fneval+1L
    }
    objp <- objlist[[vc]]$obj
    #
    if (objp<objinit) {
      lower <- vp
    } else {
      upper <- init
      vm <- max(lower, init-1L)
      vc <- as.character(vm)
      if (is.null(objlist[[vc]])) {
        objlist[[vc]] <- objfn(vm, ...)
        fneval <- fneval+1L
      }
      objm <- objlist[[vc]]$obj
      #
      if (objm<objinit) {
        upper <- vm
      } else lower <- init
    }
    #
    range <- upper - lower
    init <- floor(lower+range/2L)
    if (trace) print(c(nextlower=lower,nextupper=upper,fneval=fneval,nextinit=init))
  }
  resu <- objlist[[as.character(lower)]]
  resu$solution <- init
  resu$fneval <- fneval
  return(resu)
}

.wrap_gllim <- function(RGPpars, summstats, nbCluster, cstr=list(Sigma="i")) {
  if (length(nbCluster)>1L) {
    # Although gllim handles multiple K, 
    # it does not return enough info to extract the model with minimal AIC.
    # So each K must be fitting independently.
    objfn <- function(RGPpars,summstats, K, cstr=list(Sigma="i")) {
      cluObject <- .do_call_wrap("gllim",arglist=list(tapp=RGPpars, yapp=summstats,in_K=K, cstr=cstr), 
                                 pack="xLLiM")
      logL <- cluObject$LLf
      df <- cluObject$nbpar
      AIC <- -2*logL+2*df
      return(list(obj=AIC, cluObject=cluObject))
    }
    resu <- .wrap_minimize_over_nbClu(init=NULL, objfn=objfn, range=nbCluster, trace=FALSE, RGPpars=RGPpars, 
                      summstats=summstats, cstr=cstr)$cluObject   
  } else resu <- .do_call_wrap("gllim",arglist=list(tapp=RGPpars, yapp=summstats,in_K=nbCluster, cstr=cstr), pack="xLLiM")
  
  nclu <- length(resu$pi)
  
  RGP_COVlist <- logl_COVlist <- vector("list", nclu)
  chk <- TRUE
  for (k in seq_len(nclu)) {
    COV <- resu$Gamma[,,k] 
    # coerced to scalar if dims are 1, while drop=FALSE would keep it as array. Neither is OK , we need matrix =>
    if (is.null(dim(COV)))  dim(COV) <- c(1L,1L)
    colnames(COV) <- rownames(COV) <- rownames(RGPpars) # fittedPars
    RGP_COVlist[[k]] <- COV
    #
    COV <- resu$Sigma[,,k] 
    chk <- chk && all(diag(COV<1e-7))
    if (is.null(dim(COV)))  dim(COV) <- c(1L,1L)
    logl_COVlist[[k]] <- COV
  }
  if (chk) stop("gllim produced a degenerate resu$Sigma.")
  resu$RGP_COVlist <- RGP_COVlist
  resu$logl_COVlist <- logl_COVlist
  
  resu$clu_params$solve_t_chol_sigma_lists <- list(
    logldens=apply(resu$Sigma, 3, .solve_t_cholfn, simplify=FALSE), # of the N(y: Ak.x+bk, Sigmak) => log logLik of parms 'x' when used as function of x
    # Thus Sigma is is projected-stats space but only relevant for the 'data'
    RGPdens=lapply(RGP_COVlist, .solve_t_cholfn) # of the N(x: ck, Gammak) => marginal distrib of params in the reftable
  )
  
  class(resu) <- c("gllim",class(resu))
  resu
}

## This gives the logL of 'newRGPpars' when 'projstats'= stat.obs
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
    # Pr(k and pars) from prior Pr(k) and Pr(pars;k)
    condcluprobs[, k] <- log(object$pi[k]) + 
      .loggausspdf(t(newRGPpars), mu=object$c[, k, drop = FALSE], solve_t_chol_RGP[[k]])
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
      thr_info <- .get_thr_info(object)
      ddens <- parvaldens-thr_info$thr_dpar
      negddens <- (ddens<0)
      if (any(negddens)) {
        dlogl <- condvaldens-thr_info$reft_maxlogl
        posdlogl <- dlogl>0
        highextrapol <- negddens & posdlogl
        # make it continuous wrt to dlogl, but strongly compensating in most cases
        condvaldens[highextrapol] <- condvaldens[highextrapol]+ sqrt(dlogl[highextrapol])*ddens[highextrapol]
      }
    } else {
      stop("code missing here")
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

..simulate.gllimX <- function(object, size, parNames) {
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
  colnames(simuls) <- parNames
  return(simuls)
}

.simulate.gllimX <- function(object, 
                            seed=NULL, 
                            size=1L, # number of points for each simulation 
                            parNames, 
                            n_tables=1L,  # for # of reftables = # of bootstrap replicates
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
  if (n_tables>1L) {
    resu <- vector("list",  n_tables)
    for (it in seq_along(resu)) {
      resu[[it]] <- ..simulate.gllimX(object, size=size, parNames=parNames) 
    }
    resu
  } else ..simulate.gllimX(object, size=size, parNames=parNames) 
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

#  ____F I X M E___ hard checks for .conditional_gllimobj?
.conditional_gllimobj <- function(gllimobj, fittedPars, given, expansion=Infusion.getOption("expansion")) {  # expansion=1 to get the conditional distribution
  proportions <- gllimobj$pi
  nbCluster <- length(proportions)
  condcluprobs <- numeric(nbCluster)
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
    rhs <- try(solve(sig22,given-mean2), silent=TRUE)
    if (inherits(rhs,"try-error")) {
      sig22  <- regularize(sig22)
      rhs <- solve(sig22,given-mean2)
    }
    condmeans[clu_it,] <- MEAN[clu_it,For] + sig12 %*% rhs
    sig11 <- COV[For,For,drop=FALSE]
    RGP_COVlist[[clu_it]] <- expansion* (sig11 - sig12 %*% solve(sig22,t(sig12))) 
    condcluprobs[clu_it] <- 
      log(proportions[clu_it])+dmvnorm(t(given), # dmvnorm() tests is.vector(x) which returns FALSE if x has attributes other than names.
                                       mean2, sigma= sig22, log=TRUE)
  }
  gllimobj$c <- t(condmeans)
  gllimobj$RGP_COVlist <- RGP_COVlist
  if (is.null(dim(condcluprobs))) dim(condcluprobs) <- c(1L,length(condcluprobs))
  den <- .logsumexp(condcluprobs)
  condcluprobs <- sweep(condcluprobs, 1, den, "-")
  gllimobj$pi <- exp(condcluprobs)
  gllimobj 
}

