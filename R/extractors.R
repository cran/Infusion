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
  if (SGP) {
    .simulate.SGP(object, nsim=nsim, given=given, ...)
  } else if (inherits(object$jointdens,"dMixmod")) {
    statdens_h0 <- .conditional_Rmixmod(object$jointdens, given=given, expansion=1) 
    .simulate.MixmodResults(statdens_h0, size=nsim, drop=TRUE,
                            norm_or_t=norm_or_t) # directly in projected space
  } else if (inherits(object$jointdens,"MAF")) {
    if (is.null(dim(missing))) given <- t(given) # ugly but needed "somewhere" before reaching the Py conversion
    .simulate.MAF(object$conddens,nsim=nsim, given=given)
  } else if (inherits(object$jointdens,"dMclust")) {
    statdens_h0 <- .conditional_mclust(object$jointdens, given=given, expansion=1) 
    do.call("sim", c(statdens_h0[c("modelName", "parameters")],list(n=nsim)))[,-1,drop=FALSE] #matrix
  } else if (inherits(object$gllimobj,"gllim")) {
    .gllim.condsimul.stats(object$gllimobj, RGPpars=given, size=nsim, colTypes=object$colTypes, cbind.=FALSE)
  } 
}

