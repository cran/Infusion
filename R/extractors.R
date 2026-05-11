"get_from" <- function(object, which, ...) UseMethod("get_from") 

`get_from.default` <- function(object, which, ...) {
  return(object[[which]])
}
  
get_from.SLik <- function(object, which, raw=FALSE, force=FALSE, ...) {
  # for all versions
  if (which %in% c("obs","proj_data")) which <- "stat.obs"
  if (which %in% c("stat.obs", "Simulate", "env", "packages", 
                   "call", # this one not currently used in up-to-date workflow.
                   "control.Simulate"))  return(attr(object$logLs, which)) 
  if (raw) return(object[[which]]) # per the doc.
  ###  all (! raw) cases
  if (object$Infusion.version > "1.4.13") {
    if (which=="raw_data") { # NOT reftable.. : 'data' as opposed to 'samples'.
      return(attr(attr(object$logLs,"stat.obs"),"raw_data"))
    } else if (which=="reftable_raw") { # NOT 'data'
      return(.get_reft_raw(object))
    } else if (which=="reftable") { 
      logLs <- object$logLs
      allvars <- unlist(object$colTypes[c("inferredVars", "statNames")], use.names=FALSE)
      if ("cumul_iter" %in% colnames(logLs)) allvars <- c(allvars,"cumul_iter")
      return(logLs[,allvars]) # incl. cumul_iter
    } else if (which=="nbCluster") {
      return(.get_nbCluster_from_SLik(object))
    } else if (which=="MSL") {
      return(as.list(object[["MSL"]]))
    } else if (which=="par_RMSEs") {
      resu <- object[["par_RMSEs"]]$par_RMSEs
      if (force && is.null(resu)) {
        local_slik <- MSL(object, eval_RMSEs=TRUE, CIs=TRUE) # both needed for par_RMSEs
        local_env <- local_slik[["par_RMSEs"]]
        outer_env <- object[["par_RMSEs"]] 
        # Copy in pre-existing envir of pre-existing object: feature of get_from(), not of simple MSL() call
        for (st in ls(local_env)) assign(st,value = local_env[[st]],envir = outer_env) 
        resu <- object[["par_RMSEs"]]$par_RMSEs
      }
      return(resu)
    } else if (which=="RMSEs") {
      resu <- object[["RMSEs"]]$RMSEs
      if (force && is.null(resu)) {
        local_slik <- MSL(object, eval_RMSEs=TRUE, CIs=FALSE)
        local_env <- local_slik[["RMSEs"]]
        outer_env <- object[["RMSEs"]] 
        for (st in ls(local_env)) assign(st,value = local_env[[st]],envir = outer_env)
      }
      resu <- object[["RMSEs"]]$RMSEs
      return(resu)
    } else if (which=="C2ST") {
      misc_env <- object[["misc_env"]] # writes in this envir if it exists:
      if (is.null(resu <- misc_env$C2ST) || force) misc_env$C2ST <- resu <- .twoSampTest(object, ..., method.="C2ST")
      return(resu)
    } else if (which=="MMD") {
      misc_env <- object[["misc_env"]] # writes in this envir if it exists:
      if (is.null(resu <- misc_env$MMD) || force) misc_env$MMD <- resu <- .twoSampTest(object, ..., method.="MMD")
      return(resu)
    }
  }
  # All other cases:    # eg, reftable_raw
  return(object[[which]]) # e.g., logLs...
}

get_from.SLik_j <- get_from.SLik

logLik.SLik <- function(object, ...) structure(object$MSL$maxlogL, RMSE=unname(get_from(object,"RMSEs")[1L]))
logLik.SLik_j <- logLik.SLik


.extract_intervals <- function(object,verbose=interactive()) { ## extractor, no costly computation
  CIs <- object$CIobject$CIs
  lenCIs <- length(CIs)    
  pars <- names(CIs)
  resu <- rep(NA_real_,NROW(object$CIobject$bounds)) ## (over)size as the MSEs have no NAs
  names(resu) <- rownames(object$CIobject$bounds)
  map <- c(low=1L,up=2L)
  for (nam in rownames(object$CIobject$bounds)) {
    stterms <- regmatches(nam, regexpr("\\.", nam), invert = TRUE)[[1]] ## regmatches(<...>) splits at the first '.'
    parm <- stterms[2L]
    interval <- .get_ci_info(CIs, parm)$interval
    th <- interval[map[stterms[1L]]]
    resu[nam] <- th
  }
  if (verbose && length(resu)>0L) {
    par_headline <- "*** Interval estimates ***\n"
    cat(par_headline)
    print(resu)
    return(invisible(resu))
  } else return(resu) # may be zero-length vector if no interval info is available
}

## summary likelihood ratio (with uncertainty measures) extractor. Unfinished, in particular, need to separate residVar an to handle prior.weights
.SLR <- function(object,newdata=NULL,variance="predVar",df=NULL) {
  fittedPars <- object$colTypes$fittedPars
  if (is.null(newdata)) newdata <- unique(object$logLs[,fittedPars])
  locdata <- rbind(MSLE=object$MSL$MSLE, ## name needed for spaMM::calcNewCorrs -> newuniqueGeo
                   newdata)
  if (variance=="respVar") {
    logls <- predict(object$fit,newdata=locdata,variances=list(respVar=TRUE,cov=TRUE))
    covmat <- attr(logls,"respVar")
    # but problem with prior weights
  } else {
    logls <- predict(object$fit,newdata=locdata,variances=list(predVar=TRUE,cov=TRUE))
    covmat <- attr(logls,"predVar")
  }
  slr <- object$MSL$maxlogL -logls[]
  MSEs <- c(diag(covmat[-1,-1,drop=FALSE])+covmat[1,1]-2*covmat[1,-1])
  if ( any(MSEs<0) ) {
    message("Inaccurate MSE computation, presumably from nearly-singular covariance matrices.")
    MSEs[MSEs<0] <- NA ## quick patch
  }
  RMSEs <- sqrt(MSEs)
  resu <- rbind(LRstat=2*slr[-1],RMSE=RMSEs)
  if (!is.null(df)) {
    resu <- rbind(resu,p_RMSE=RMSEs*dchisq(resu["LRstat",],df=df))
  }
  return(resu)
}

#.SLR(slik)

"summLik" <- function(object, parm, data, ...) UseMethod("summLik") ## makes it easy to develop new inference methods

summLik.SLik_j <- function(object, parm, data=t(get_from(object,"proj_data")), log=TRUE, 
                           which="safe", constr_tuning = Inf, newMSL=FALSE, ...) {
  fittedPars <- object$colTypes$fittedPars
  is1Dprof <- is.null(dim(parm)) && # profile performed for single parm vector only. 
    # By contrast a 'parm' table must contain all fittedPars. See comment on recursive call below.
    length(profiledOutpars <- setdiff(fittedPars, names(parm)))
  if (is1Dprof) {
    # prof_init_pars <- object$MSL$MSLE[profiledOutpars] 
    prof_lower <- object$lower[profiledOutpars]
    prof_upper <- object$upper[profiledOutpars]
    
    ## based on .ecdf_2lr:
    # Maximization of likelihood wrt profiledOutpars (for the newobs):
    # it should be like running profile(object, value=BGP[h0_pars]) except that new data are used hence not the lik surf of the 'object'
    if ( ! missing(data)) {
      if (newMSL) { # note default is FALSE. 
        # This block evaluates opt_mlogL_Dsim used to eval attributes optionally added added to the profile. 
        init_h1 <- .safe_init(object=object, newobs = data, plower=object$lower,pupper=object$upper)
        # Next result is "ML" for the newdata (not prof lik for tested hypo) => there's an optim for all params in object$lower.
        opt_mlogL_Dsim <- .optim_mlogL_newobs(object, newobs=data,
                                              init=unlist(init_h1), 
                                              lower=object$lower, upper=object$upper,
                                              which=which)
        init_h0 <- .safe_init(object, given=parm, plower=prof_lower, pupper=prof_upper,
                              newobs=data, MSLE=opt_mlogL_Dsim$solution, # controls call of .get_init_from_hull()
                              profiledNames=profiledOutpars) 
      } else init_h0 <- .safe_init(object, given=parm, plower=prof_lower, pupper=prof_upper,
                                   newobs=data, MSLE=NULL, # explicit NULL important here (no call of .get_init_from_hull())
                                   profiledNames=profiledOutpars) # accounts for constr_crits...
    } else { # no new data
      init_h0 <- .safe_init(object, given=parm, plower=prof_lower, pupper=prof_upper,
                            profiledNames=profiledOutpars) # calling .get_init_from_hull()
    }
    opt_prof_mlogL_Dsim <- .optim_mlogL_newobs(object, 
                                               givenpars=parm, # would be useful only for non-default (object$MSL$MSLE) values
                                               newobs=data, 
                                               init=unlist(init_h0)[profiledOutpars], 
                                               lower=prof_lower, upper=prof_upper,
                                               which=which)
    logL_h0 <- - opt_prof_mlogL_Dsim$objective
    profpt <- c(opt_prof_mlogL_Dsim$solution, parm)[object$colTypes$fittedPars]
    if (exists("opt_mlogL_Dsim")) {
      attr(logL_h0,"LRstat") <- - opt_mlogL_Dsim$objective - logL_h0
      attr(logL_h0,"newobs_MSL") <- list(MSLE=opt_mlogL_Dsim[["solution"]],
                                         maxlogL=-opt_mlogL_Dsim[["objective"]])
    }
    attr(logL_h0,"profpt") <- profpt
    logL_h0
  } else { # not 1Dprof
    ## Note recursive call from
    ## summLik.SLik_j(., incomplete parm) -> .optim_mlogL_newobs -> .safe_init -> summLik.SLik_j(., parm =<table of full-par inits> )
    predict.SLik_j(object, 
                   newdata=parm, ## requests full new fittedPars values! -- may have multiple rows
                   log=log, 
                   which=which, # may be still preferable, to "safe", for drawing new points
                   tstat= data, # 1-row matrix...
                   constr_tuning=constr_tuning,
                   ...)
    
  }
}

# has SLikp and SLik method

`summLik.default` <- function(object, parm, data=t(get_from(object,"proj_data")), ...) predict(object, newdata=parm, tstat=data, ...)

# private bc no use currently.
.hessian.summLik <- function(fitobject,
                            x=c(fitobject$MSL$MSLE, get_from(fitobject,"proj_data")), 
                            which="safe", ...) {
  colTypes <- fitobject$colTypes
  fittedPars <- colTypes$fittedPars
  statNames <- colTypes$statNames
  locfn <- function(z) {
    summLik.SLik_j(fitobject, parm=z[fittedPars], data=t(z[statNames]), which=which, ... )
  }
  hess <- hessian(locfn, x=x) # again pracma::hessian may give better results.
  colnames(hess) <- rownames(hess) <- names(x)
  hess
}

.simulate.SGP <- function(object, nsim=1L, given,
                          Simulate=get_from(object,"Simulate"),
                          control.Simulate=get_from(object,"control.Simulate"),
                          ...) {
  if (! (is.function(Simulate) || is.character(Simulate))) 
    stop("Simulate() function not available from SLik object.")
  .wrap_Simulate(Simulate,  
                 parsTable=data.frame(t(given))[rep(1,nsim),], 
                 control.Simulate=control.Simulate,
                 reftable=NULL, ...)
}

simulate.SLik_j <- function(object, nsim = 1L, 
                            seed = NULL, # seed ignored in all cases: rethink ____F I X M E____
                            given=object$MSL$MSLE, 
                            norm_or_t=.wrap_rmvnorm, 
                            SGP=FALSE, ...) {
  if (SGP) { # sample-generating process
    .simulate.SGP(object, nsim=nsim, given=given, ...)
  } else if (inherits(object$jointdens,"dMixmod")) {
    statdens_h0 <- .conditional_Rmixmod(object$jointdens, given=given, expansion=1) 
    .simulate.MixmodResults(statdens_h0, size=nsim, drop=TRUE,
                            norm_or_t=norm_or_t) # directly in projected space
  } else if (inherits(object$jointdens,"MAF")) {
    # If using="c.mafR", object$conddens is a MAF, else it may be NULL; 
    if (inherits(object$conddens,"MAF")) {
      if (is.null(dim(given))) given <- t(given) # ugly but needed "somewhere" before reaching the Py conversion
      .simulate.MAF(object$conddens,nsim=nsim, given=given) 
    } else if (is.null(object$conddens)) {
      # If using="MAFmix", there is no $conddens
      # MAFmix provides two MAFs, $jointdens and $pardens, and a gaussian mixture, $MGMjoindens.
      conddens <- .get_conditional_GMM(object,given=given)
    } else stop("Case not yet considered.")
  } else if (inherits(object$jointdens,"dMclust")) {
    statdens_h0 <- .conditional_dMclust(object$jointdens, given=given, expansion=1, using=object$using) 
    if (object$using=="mclust") {
      do.call("sim", c(statdens_h0[c("modelName", "parameters")],list(n=nsim)))[,-1,drop=FALSE] #matrix
    } else .simVVV_rmvnorm_or_t(parameters=statdens_h0$parameters,n=nsim, use_rmvt=FALSE)[,-1,drop=FALSE] #matrix
  } else if (inherits(object$gllimobj,"gllim")) {
    .gllim.condsimul.stats(object$gllimobj, RGPpars=given, size=nsim, colTypes=object$colTypes, cbind.=FALSE)
  } 
}

.mmd_biased <- function(XX, YY, XY, n1=nrow(XX), n2=nrow(YY)){
  sum(XX)/(n1^2) + sum(YY)/(n2^2) - (2/(n1*n2))*sum(XY)
}

.mmd_unbiased <- function(XX, YY, XY, n1=nrow(XX), n2=nrow(YY)){
  a <- (sum(XX)-sum(diag(XX)))/(n1*(n1-1L))
  b <- (sum(YY)-sum(diag(YY)))/(n2*(n2-1L))
  c <- (2/(n1*n2))*sum(XY)
  a+b+c
}

..MMDtest <- function(mat, n1, n2, method=c("b","u"), nperm) {
  if (n1<2L) stop("n1 sould be >1 !")
  if (n2<2L) stop("n2 sould be >1 !")
  method <- match.arg(tolower(method), c("b","u"))
  id1 <- seq(n1)
  ntot <- n1+n2
  id2 <- (n1+1L):ntot
  statfn <- switch(method,
                   "b" = .mmd_biased,
                   "u" = .mmd_unbiased)
  
  stat0 <- statfn(mat[id1,id1], mat[id2,id2], mat[id1,id2],n1=n1,n2=n2)
  
  statvec <- numeric(nperm)
  for (i in seq_len(nperm)){
    perm <- sample(ntot)
    pid1 <- perm[id1]
    pid2 <- perm[id2]
    statvec[i] <- statfn(mat[pid1,pid1], mat[pid2,pid2], mat[pid1,pid2],n1=n1,n2=n2)
  }
  pvalue <- (sum(statvec>=stat0)+1L)/(nperm+1L)
  list(p_value=pvalue, n1=n2, n2=n2, method=method, nperm=nperm)
}

# writes in misc_env as side effect
.twoSampTest <- function(object, nsim=1000L, nperm=999L,  
                         method.="C2ST", resimulate=FALSE, 
                         Simulate=get_from(object,"Simulate"), 
                         method="b", # for MMD
                         check=FALSE
                         ) {
  statNames <- object$colTypes$statNames
  fittedPars <- object$colTypes$fittedPars
  nsim <- min(nrow(object$logLs), nsim)
  X <- tail(object$logLs, nsim)
  parsTable <- X[,fittedPars, drop=FALSE]
  if (check) { # check unif distrib of p-values... 
    # actually oob in C2ST then "under classifies": deficiency of low p-values
    # so when oob is used in the real compar between Simulate() and simulate(), 
    # and the two distribs differ, one may hae a sigmoid ecdf.
    X <- vector("list", nsim)
    for (it in seq_len(nsim))  X[[it]] <- simulate(object,1, given=unlist(parsTable[it,]))
    X <- data.frame(do.call(rbind,X))
    X <- cbind(parsTable, X)
  } else if (resimulate) {
    trypoints <- cbind(parsTable, object$colTypes$fixedPars) ## add fixedPars for simulation
    trypoints <- trypoints[,object$colTypes$allPars,drop=FALSE] ## column reordering
    X <- add_reftable(, Simulate = Simulate, parsTable = trypoints)
    if (length(object$projectors)) X <- project(X, projectors=object$projectors)
  } 
  dgpS <- X[,c(fittedPars,statNames), drop=FALSE]
  #
  fitS <- vector("list", nsim)
  for (it in seq_len(nsim))  fitS[[it]] <- simulate(object,1, given=unlist(parsTable[it,]))
  fitS <- data.frame(do.call(rbind,fitS))
  fitS <- cbind(parsTable, fitS)
  #
  if (method.=="C2ST") {
    dgpS$.class <- "mgm"
    fitS$.class <- "sample"
    train <- rbind(fitS,dgpS)
    train$.class <- as.factor(train$.class) 
    form <- as.formula(paste(".class ~ ",paste(statNames,collapse="+")))
    rf <- ranger(formula = form, data = train)
    conf <- rf$confusion.matrix
    object$misc_env$C2ST <- 
      c2st <- list(test=binom.test(sum(diag(conf)),n=2*nsim,alternative = "greater"))
    c2st
  } else if (method.=="MMD") {
    mmdsample <- rbind(fitS,dgpS)
    distmat <- as.matrix(proxy::dist(mmdsample))
    mat <- exp(-distmat^2)   
    object$misc_env$MMD <- 
      mmdtest <- ..MMDtest(mat, n1=nsim, n2=nsim, method=method, nperm=nperm)
    mmdtest
  } else stop("Unknown 'method'.")
}

C2ST <- function(
    object, nsim, 
    resimulate=FALSE,
    Simulate=get_from(object,"Simulate"), 
    ...) .twoSampTest(
      object=object, nsim=nsim, resimulate=resimulate, Simulate=Simulate,
      method.="C2ST", ...)