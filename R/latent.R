declare_latent <- function(reftable, latentVars) {
  attr(reftable,"latentVars") <- latentVars
  attr(reftable,"LOWER") <- c(attr(reftable,"LOWER"), 
                              setNames(rep(-Inf, length(latentVars)), latentVars ))
  attr(reftable,"UPPER") <- c(attr(reftable,"UPPER"), 
                              setNames(rep(Inf, length(latentVars)), latentVars ))
  reftable
}

pplatent <- function(object, type="mean",
                     newDP=NULL,
                     sumstats= t(get_from(object,"stat.obs")), # default= 1-row matrix...
                     pars=t(object$MSL$MSLE), # default= 1-row matrix...
                     size=1000L,
                     ...) {
  if (is.null(dim(sumstats))) sumstats <- t(sumstats) # matrix needed below
  if (is.null(newDP)) newDP <- as.matrix(cbind(as.data.frame(pars),as.data.frame(sumstats))) # ugly but it works... 
  # 'given' must be vector not data frame
  nr <- nrow(newDP)
  latentVars <- object$colTypes$latentVars
  pred <- matrix(NA, nrow=nr, ncol=length(latentVars))
  
  if (inherits(object$jointdens,"dMixmod")) {
    for (ii in seq_len(nr)) {
      latentdens <- .conditional_Rmixmod(object$completedens, given=newDP[ii,])
      if (type=="median") {
        pseudoVs <- .simulate.MixmodResults(latentdens, size=size, drop=TRUE,
                                            norm_or_t=.wrap_rmvnorm)
        pred[ii,] <- matrixStats::colMedians(pseudoVs)
      } else pred[ii,] <- matrixStats::colSums2(latentdens@parameters@proportions * latentdens@parameters@mean)
    }
  } else if (inherits(object$jointdens,"dMclust")) {
    for (ii in seq_len(nr)) {
      latentdens <- .conditional_mclust(object$completedens, given=newDP[ii,])
      if (type=="median") {
        pseudoVs <- do.call("sim", c(latentdens[c("modelName", "parameters")],list(n=size)))[,-1,drop=FALSE] #matrix
        pred[ii,] <- matrixStats::colMedians(pseudoVs)
      } else pred[ii,] <- matrixStats::colSums2(latentdens$parameters$pro * t(latentdens$parameters$mean))
    }
  } else if (inherits(object$gllimobj,"gllim")) { # speculative, never tested
    fittedPars <- object$colTypes$fittedPars 
    for (ii in seq_len(nr)) {
      latentdens <- .conditional_gllimobj(object$gllimobj, fittedPars=fittedPars, 
                                          given=newDP[ii,])
      if (type=="median") {
        pseudoVs <- .simulate.gllimX(latentdens, size=size, parNames=object$colTypes$latentVars)
        pred[ii,] <- matrixStats::colMedians(pseudoVs)
      } else pred[ii,] <- matrixStats::colSums2(latentdens$pi * t(latentdens$c)) # t() part of the speculation
    }
  }
  
  pred <- drop(pred)
  names(pred) <- latentVars
  pred
}

.platent_1DP <- function(q, # expects vector, 1-row matrix should work [vectors of latent vars for one sample of sample-GP] 
                         latentdens,
                         proportions,
                         means,
                         sds # sds is here a *list*
                         ) {
  nbClu <- length(proportions)
  CDF <- matrix(NA, ncol=length(q), nrow=nbClu)
  for (clu_it in seq_len(nbClu)) {
    CDF[clu_it, ] <- proportions[clu_it] * pnorm(q,means[clu_it,], sd=sds[[clu_it]])
  }
  CDF <- matrixStats::colSums2(CDF)
  CDF # vector
}

# and the \code{platent} and \code{qlatent} are the distribution function  
# of the quantile function for the latent variable given the fitted parameters 
# and data. Using \code{qlatent} to obtain prediction intervals 
# is the naive plug-in method that ignores uncertainty in parameter estimates. 

.platent <- function(q, # expects vector, 1-row matrix should work [vectors of latent vars for one sample of sample-GP] 
                     object, newDP=NULL,
                    sumstats= t(get_from(object,"stat.obs")), # default= 1-row matrix...
                    pars=t(object$MSL$MSLE), # default= 1-row matrix...
                    ...) {
  if (is.null(dim(sumstats))) sumstats <- t(sumstats) # matrix needed below
  if (is.null(newDP)) newDP <- as.matrix(cbind(as.data.frame(pars),as.data.frame(sumstats))) # ugly but it works... 
  # 'given' must be vector not data frame
  nr <- nrow(newDP)
  CDF <- matrix(NA, nrow=nr, ncol=length(q))
  completedens <- object$completedens
  if (inherits(completedens,"dMixmod")) {
    for (ii in seq_len(nr)) {
      latentdens <- .conditional_Rmixmod(completedens, given=newDP[ii,])
      parameters <- latentdens@parameters
      sds <- lapply(parameters@variance, function(v) sqrt(diag(v))) # only considering marginal distribus of each latent var
      CDF[ii,] <- .platent_1DP(q, latentdens,
                               proportions=parameters@proportions,
                               means=parameters@mean,
                               sds= sds ) # sds is here list
    }
  } else if (inherits(completedens,"dMclust")) {
    for (ii in seq_len(nr)) {
      latentdens <- .conditional_mclust(completedens, given=newDP[ii,])
      parameters <- latentdens$parameters
      sds <- apply(parameters$variance$sigma,3L, function(v) sqrt(diag(v)), simplify=FALSE)
      CDF[ii,] <- .platent_1DP(q, latentdens,
                               proportions=parameters$pro,
                               means=t(parameters$mean),
                               sds=sds  ) # idem
    }
  }
  
  CDF
}

.qlatent_1DP <- function(p, latentdens, 
                         s_latentv, 
                         means, 
                         sds, # reduced to a vector...
                         proportions) {
  # s_latentv are latent values in the reftable, used to sample (in an uncontrolled way, but widely)
  # the range of the distribution of the latent variable, as a way of avoid numerical root finding. 
  # The precise distribution of s_latentv should have little impact on the result provided sampling is dense.
  CDF <- 0
  for(k in seq_along(means)) {
    CDF <- CDF + proportions[k] * pnorm(s_latentv,means[k], sds[k]) 
    # => a value in (0,1) _for each of the many elements of s_latentv_ depending on its value relative relative to means
    # 's'_latent values are sorted in increasing order so CDF values are sorted in increasing order too
  }
  quant <- numeric(length(p))
  for (p_it in seq_along(p)) {
    # On (CDF) values sorted in increasing order, use which.max(.>.) or which.min(.<.) ...  :
    # Counterintuitively, it returns the first of the sorted (CDF) values that matches the condition.
    # which.max is then the fastest method to find first value that matches the condition
    whichsup <- which.max(CDF>p[p_it]) 
                                       
    minsup <- min(whichsup)
    if(minsup==Inf) { # unexpected, but not caused by forgetting to sort s_latentv   # no CDF value higher than p
      quant[p_it] <- Inf 
    } else if(minsup==1L) { 
      if (any( ! whichsup)) stop(".qlatent_1DP() called on unsorted 's_latentv'!")
      quant[p_it] <- -Inf # all CDF values higher than p
    } else quant[p_it] <- sum(s_latentv[minsup+c(-1L,0L)])/2 # This is where root finding would be needed.
  }
  quant
}

.qlatent_1latent_Rmixmod <- function(
    s_latentv, locdens, newDPs, 
    p, # percentiles 'Gw', corrected relative to nominal ones, to approximate nominal (intended) coverage of the result 
    levels # The intended coverages, used here only as colnames for return object.
) {
  s_latentv <- sort(s_latentv)
  nr <- nrow(newDPs)
  quant <- matrix(NA, nrow=nr, ncol=length(p))
  for (ii in seq_len(nr)) {
    latentdens <- .conditional_Rmixmod(locdens, given=newDPs[ii,]) # must be 1D
    # In default calls, newDPs is simply SMLE+prod_data so latentdens is only considered in such 'best' conditions 
    parameters <- latentdens@parameters
    # Get the quantiles of the latent var for the given percentile,
    #    'p', which has been corrected relative to the nominal coverage probability.
    quant[ii,] <- .qlatent_1DP(
      p=p, latentdens, 
      s_latentv=s_latentv, # again, only used as a dense coverage of the range, not fortheir precise distribution.
      means=parameters@mean,
      sds=sqrt(unlist(parameters@variance,use.names = FALSE)), # (single vector as each matrix is a single variance)  
      proportions=parameters@proportions
    )
  }
  colnames(quant) <- levels
  quant
}

.qlatent_1latent_mclust<- function(s_latentv, locdens, newDPs, p, 
                             levels # colnames for return object
) {
  s_latentv <- sort(s_latentv)
  nr <- nrow(newDPs)
  quant <- matrix(NA, nrow=nr, ncol=length(p))
  for (ii in seq_len(nr)) {
    latentdens <- .conditional_mclust(locdens, given=newDPs[ii,])
    parameters <- latentdens$parameters 
    quant[ii,] <- .qlatent_1DP(p=p, latentdens, s_latentv=s_latentv,
                               means=t(parameters$mean),
                               sds=apply(parameters$variance$sigma,3,sqrt), # single vector...
                               proportions=parameters$pro)
  }
  colnames(quant) <- levels
  quant
}


.qlatent <- function(p, object, focalDPs=NULL,
                    sumstats= t(get_from(object,"stat.obs")), # default= 1-row matrix...
                    pars=t(object$MSL$MSLE), # default= 1-row matrix...
                    levels # colnames for return object
                    ) {
  if (is.null(dim(sumstats))) sumstats <- t(sumstats) # matrix needed below
  statNames <- names(sumstats)
  if (is.null(focalDPs)) {
    # ugly but it works (given' must be vector not data frame):
    focalDPs <- as.matrix(cbind(as.data.frame(pars),as.data.frame(sumstats))) # 1-row for default pars 
  }  
  latentVars <- object$colTypes$latentVars
  s_latentv <- object$logLs[,latentVars, drop=FALSE]
  # s_latentv is only expected to be a dense exploration of the range of latent values,
  # this property being used in .qlatent_1latent_Rmixmod() to avoid numerical root finding by a simple approx.
  # The distrib of s_latentv (as determined by the SGPs i nthe reftable) should then not matter.
  # .qlatent_1latent_Rmixmod() -> .qlatent_1DP() expect s_latentv to be sorted!
  # But sorting can only be performed when a single column (single latentVar) has been selected.
  completedens <- object$completedens
  resu <- vector("list",length(latentVars))
  names(resu) <- latentVars
  if (inherits(completedens,"dMixmod")) {
    varNames <- object$completedens@varNames
    for (st in latentVars) {
      # Perhaps not essential step, slightly complicating the understanding of the code:
      locdens <- .marginalize_Rmixmod(object$completedens, colNames=varNames ,
                                      For=setdiff(varNames, setdiff(latentVars,st)))
      # locdens is NOT YET the conditional distrib of the latent var given the focalDP.
      # It retains all dimensions in 'For', thus only latentVars not currently inferred are removed.
      # This .marginalize...() was conceived as an optimization 
      # of later .conditional...()'s for hypothetical different rows of focalDP. 
      resu[[st]] <- .qlatent_1latent_Rmixmod(s_latentv=s_latentv[,st], locdens=locdens, 
                                             newDPs=focalDPs, p=p, levels=levels)
    }
  } else if (inherits(completedens,"dMclust")) {
    varNames <- rownames(completedens$parameters$mean)
    for (st in latentVars) {
      locdens <- .marginalize_mclust(object$completedens, colNames=varNames ,
                                      For=setdiff(varNames, setdiff(latentVars,st)))
      resu[[st]] <- .qlatent_1latent_mclust(s_latentv=s_latentv[,st], locdens=locdens, 
                                            newDPs=focalDPs, p=p, levels=levels)
    }
  }
  resu
}

.dlatent_1DP <- function(x=x, latentdens) {
  parameters <- latentdens@parameters
  sds <- sqrt(unlist(parameters@variance,use.names = FALSE))
  sum(parameters@proportions * dnorm(x,parameters@mean, sd=sds))
}

# dlatent not used -> latest code not checked and not extended to dMclust objects. 
.dlatent <- function(x, object, newDP=NULL,
                     sumstats= t(get_from(object,"stat.obs")), # default= 1-row matrix...
                     pars=t(object$MSL$MSLE), # default= 1-row matrix...
                     ...) {
  if (is.null(dim(sumstats))) sumstats <- t(sumstats) # matrix needed below
  if (is.null(newDP)) newDP <- as.matrix(cbind(as.data.frame(pars),as.data.frame(sumstats))) # ugly but it works... 
  # 'given' must be vector not data frame
  nr <- nrow(newDP)
  PDF <- numeric(nr)
  for (ii in seq_len(nr)) {
    latentdens <- .conditional_Rmixmod(object$completedens, given=newDP[ii,]) # all functions assume a single latentVar
    # otherwise one has to condition on 'other' latentVars in order to get a one dimensional distrib
    parameters <- latentdens@parameters
    PDF[ii] <- .dlatent_1DP(x=x, latentdens)
  }
  PDF
}

# Implementation of LawlessF05's method
latint <- function(object, nsim=199L, levels=c(0.025,0.975), 
                   sumstats= t(get_from(object,"stat.obs")), # default= 1-row matrix...
                   Simulate,
                   control.Simulate=get_from(object,"control.Simulate"),
                   bootSamples=NULL,
                   ...) {
  latentVars <- object$colTypes$latentVars
  statNames <- object$colTypes$statNames
  
  # 1) simulate bootstrap replicates of (data and latentvar) D* and O*
  if (missing(Simulate)) {
    Simulate <- get_from(object,"Simulate")
    if (is.null(Simulate)) warning(paste("Simulate() function not available from SLik object.\n",
                                   "Samples will be simulated using the gaussian mixture model,\n",
                                   "which appears to give less reliable intervals."), immediate. = TRUE)
    # The user an easily suppress the warning by providing an explicit NULL
  }
  if (is.function(Simulate) || is.character(Simulate)) { # RECOMMENDED DEFAULT
    #newsimuls <- t(replicate(nsim, Simulate(parvec=object$MSL$MSLE))) # parvec syntax not universal ?
    #newsimuls <- t(replicate(nsim, do.call(Simulate, as.list(object$MSL$MSLE)))) 
    if (is.null(bootSamples)) bootSamples <- .wrap_Simulate(
      Simulate,  
      parsTable=data.frame(t(object$MSL$MSLE))[rep(1,nsim),], 
      control.Simulate=control.Simulate,
      reftable=NULL, ...)
    latentsimuls <- bootSamples[,latentVars]
    if (is.null(object$projectors)) {
      ranprojSamples <- bootSamples
    } else ranprojSamples <- .project_reftable_raw(bootSamples, projectors = object$projectors, ext_projdata=object$projdata)
  } else {    # code tidying not tested since this option is not recommended
    # directly in projected space en keeping the latent variable
    ranprojSamples <- simulate(object, nsim=nsim) # n replicates sample-GP
    latentsimuls <- ranprojSamples[,latentVars] # nrow = n replicates sample-GP
  }
  newprojStats <- ranprojSamples[,statNames, drop=FALSE] # nrow = n replicates sample-GP
  
  latentsimuls <- as.matrix(latentsimuls) # bc  simfun ->.platent -> .platent_1DP -> pnorm does not handle data.frame
  
  simfun <- function(it) { # other arguments are in the environment of the function, which seems sufficient.
    # 2) refit the model to data D* to obtain theta(D*)
    # 3) evaluate p* = CDF p...(new latent var, theta(D*)) 
    newobs <- newprojStats[it,, drop=FALSE] # data frame for one sample from sample-GP
    # Maximization of likelihood wrt full param (for the simulated data):
    #init <- BGP
    opt_mlogL_Dsim <- .optim_mlogL_newobs(object, newobs=newobs,
                                          init=NULL, which="safe",
                                          lower=object$lower, upper=object$upper)
    .platent(q = latentsimuls[it,], # vector of latent vars for one replicate of sample-GP
            sumstats = newobs, # data frame of projections (for params and latent vars) for one replicate of sample-GP
            object = object, pars=t(opt_mlogL_Dsim$solution))
  }
  iterator <- seq_len(nsim)
  dim(iterator) <- c(1L,nsim)
  object <- .shrink(object, ...) 
  pboo <- spaMM::dopar(newresp = iterator, fn = simfun, fit_env=list(), ... )
  
  # 4) evaluate G(w) = bootstrap average of I(p*<w) for different w (says 0.025 and 0.975)
  Gw <- sapply(levels, function(w) sum(pboo<=w)/length(pboo)) # 0.02120212 0.96019602 which looks consistent with the plot
  
  # 5) and presumably invert this using q(., theta(D)) to obtain interval on response scale.
  .qlatent(p=Gw, sumstats=sumstats, object=object, levels=levels)
}
