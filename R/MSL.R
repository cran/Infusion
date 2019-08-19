.RMSEwrapper <- function(object, CIpoints=object$CIobject$bounds, useCI=TRUE, verbose=interactive()) {
  if( useCI && ! is.null(CIpoints) ) {
    locdata <- data.frame(rbind(MSLE=object$MSL$MSLE,CIpoints))
    etas <- predict(object, newdata=locdata,variances=list(linPred=TRUE,dispVar=TRUE,cov=TRUE))
    covmat <- attr(etas,"predVar")
    #object$CIobject$cov <- covmat
    MSEs <- c(MSL=covmat[1,1],diag(covmat[-1,-1,drop=FALSE])+covmat[1,1]-2*covmat[1,-1])
  } else {
    etas <- predict(object,newdata=object$MSL$MSLE,variances=list(linPred=TRUE,dispVar=TRUE))
    MSEs <- structure(attr(etas,"predVar")[1L],names="MSL") ## of predict= of logL
  }
  if ( any(MSEs<0) ) {
    message("Inaccurate MSE computation, presumably from nearly-singular covariance matrices.")
    # in particular in infer_surface.tailp 
    # infinite lambda => MSE < 0 => nvec[1]<0 => size<0 in rvolTriangulation crashes 
    MSEs[MSEs<0] <- NA ## quick patch
  }  
  if (inherits(object,"SLikp")) {
    dmudeta <- object$fit$family$mu.eta(as.numeric(etas)) ## also removes attributes, quite useful
    MSEs <- MSEs * (dmudeta^2) 
  }
  RMSEs <- sqrt(MSEs)
  if (length(MSEs)>1L) {
    if (inherits(object,"SLik")) {
      headline <-  "* RMSEs of summary log-L maximum and of its ratio at CI bounds: *\n"
    } else if (inherits(object,"SLikp")) {
      headline <-  "* RMSEs of MaxSumm-tail p and of MSQ for CIs: *\n"
    }
  } else {
    if (inherits(object,"SLik")) {
      headline <-  "* RMSEs of summary log-L maximum: *\n"
    } else if (inherits(object,"SLikp")) {
      headline <-  "* RMSE of MaxSumm-tail p: *\n"
    }      
  }
  if(verbose) {
    cat(headline)
    print(RMSEs)
  }
  return(RMSEs)
}

# creates <slik>$par_RMSEs 
.par_RMSEwrapper <- function(object,verbose=interactive()) {
  RMSEs <- object$RMSEs
  if( (lenRMSEs <- length(RMSEs))>1L)  {
    CIs <- object$CIobject$CIs
    pars <- names(CIs)
    par_RMSEs <- rep(NA,lenRMSEs-1L) ## (over)size as the MSEs have no NAs
    par_ests <- rep(NA,lenRMSEs-1L) ## (over)size as the MSEs have no NAs
    names(par_RMSEs) <- names(par_ests) <- names(RMSEs)[-1L]
    map <- c(low=1L,up=2L)
    for (nam in names(par_RMSEs)) {
      stterms <- regmatches(nam, regexpr("\\.", nam), invert = TRUE)[[1]] ## regmatches(<...>) splits at the first '.'
      parm <- stterms[2L]
      interval <- CIs[[parm]]$interval
      th <- interval[map[stterms[1L]]]
      par_ests[nam] <- th
      if (! is.na(RMSEs[nam])) {
        objfn <- function(value) {
          names(value) <- parm
          profile(object, value=value)
        } 
        dlogLdth <- grad(objfn,x=th)
        par_RMSEs[nam] <- RMSEs[nam]/abs(dlogLdth)
      } else par_RMSEs[nam] <- NA
    }
    resu <- rbind(par=par_ests,par_RMSE=par_RMSEs,LR_RMSE=(RMSEs)[-1L]) 
    if (verbose && length(par_RMSEs)>0L) {
      par_headline <- "*** Interval estimates and RMSEs ***\n"
      cat(par_headline)
      print(resu)
      return(invisible(resu))
    } else return(resu)
  } else return(NULL)
}

.par_wrapper <- function(object,verbose=interactive()) {
  CIs <- object$CIobject$CIs
  lenCIs <- length(CIs)    
  pars <- names(CIs)
  resu <- rep(NA_real_,NROW(object$CIobject$bounds)) ## (over)size as the MSEs have no NAs
  names(resu) <- rownames(object$CIobject$bounds)
  map <- c(low=1L,up=2L)
  for (nam in rownames(object$CIobject$bounds)) {
    stterms <- regmatches(nam, regexpr("\\.", nam), invert = TRUE)[[1]] ## regmatches(<...>) splits at the first '.'
    parm <- stterms[2L]
    interval <- CIs[[parm]]$interval
    th <- interval[map[stterms[1L]]]
    resu[nam] <- th
  }
  if (verbose && length(resu)>0L) {
    par_headline <- "*** Interval estimates ***\n"
    cat(par_headline)
    print(resu)
    return(invisible(resu))
  } else return(resu)
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
    logls <- predict(object$fit,newdata=locdata,variances=list(linPred=TRUE,dispVar=TRUE,cov=TRUE))
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

# both SLik and SLikp, with different methods used in -> allCIs -> confint
MSL <- function (object,CIs=TRUE,level=0.95, verbose=interactive(),
                 eval_RMSEs=TRUE, #inherits(object,"SLik"), 
                 ...) { ##
  # Maximization ## revised 13/07/2016
  fittedPars <- object$colTypes$fittedPars
  if (inherits(object,"SLik")) {
    vertices <- object$fit$data[,fittedPars,drop=FALSE]
    lowup <- apply(vertices,2L,range)
    lower <- lowup[1L,]
    upper <- lowup[2L,]
  } else {
    if (verbose) cat(paste(nrow(object$logLs),"simulated samples;\n"))
    lower <- object$lower
    upper <- object$upper
  }
  # (a constrOptim woul be even better)
  stat.obs <- object$stat.obs
  # coeffs <- object$coeffs
  parscale <- (upper-lower) 
  if ( is.null(init <- object$optrEDF$par)) { ## ...if inherits(object,"SLik_j")
    # SLik object:
    init <- unlist(object$obspred[which.max(object$obspred[,attr(object$obspred,"fittedName")]),fittedPars]) 
    init <- init*0.999+ colMeans(vertices)*0.001
  }
  method <- "L-BFGS-B" ## works also in 1Dand does not ignore the init value (while Brent would)
  time1 <- Sys.time()
  msl <- optim(init, function(v) {as.numeric(predict(object,newdata=v))},
               ## as numeric because otherwise in 1D, optim -> minimize -> returns a max 
               ##   of same type as predict(object$fit,newdata=v), i.e. matrix... 
               ## FR->FR with spaMM>1.7.12, predict(object,newdata=v)[] should be OK  
               lower=lower,upper=upper,control=list(fnscale=-1,parscale=parscale),method=method)
  optim_time <- round(as.numeric(difftime(Sys.time(), time1, units = "secs")), 1) ## spaMM:::.timerraw(time1)
  if (inherits(object,"SLik") && length(fittedPars)>1L) {
    locchull <- resetCHull(vertices, formats=c("constraints"))
    if( ! (isPointInCHull(msl$par, constraints=locchull[c("a", "b")]))) { ## if simple optim result not in convex hull
      ui <- -locchull$a
      ci <- -locchull$b
      objectivefn <- function(v) {as.numeric(predict(object,newdata=v))}
      objectivefn.grad <- function(x) {grad(func=objectivefn, x=x)} ## no need to specify fixedlist here
      msl <- constrOptim(theta=init,f=objectivefn,grad=objectivefn.grad, ui=ui, ci=ci , 
                         mu=1e-08, ## a low mu appear important
                         method = "BFGS",control=list(fnscale=-1,parscale=parscale))
      ## then we will optimize in the convex hull using constroptim(R). But we need a good starting point
    }
  }
  if (verbose) {
    if (inherits(object,"SLik")) {
      cat("* Summary ML: *\n")
      print(c(msl$par,"logL"=msl$value))
    } else if (inherits(object,"SLik_j")) {
      cat("* Summary ML: *\n")
      print(c(msl$par,"logL"=msl$value))
    } else if (inherits(object,"SLikp")) {
      cat("* Summary MQ: *\n")
      print(c(msl$par,"tailp"=msl$value))
    }
  }
  object$MSL$MSLE <- msl$par
  if (length(lower)==1) names(object$MSL$MSLE) <- names(lower)
  object$MSL$maxlogL <- msl$value
  if (inherits(object$fit,"HLfit")) object$MSL$predVar <- attr(predict(object,newdata=msl$par,variances=list(linPred=TRUE,dispVar=TRUE)),"predVar")
  # CIs
  if(CIs) {
    locverbose <- (verbose && ! inherits(object$fit,"HLfit")## for HLfit object, printing is later
                   && optim_time>3) # guess when it is useful to be verbose from the time to find the maximum 
    if (locverbose) {
      prevmsglength <- .overcat("Computing confidence intervals...\n", 0L) 
    } else {
      prevmsglength <- 0L
    }
    object$CIobject <- .allCIs(object,verbose=locverbose, level=level)  ## may be NULL
  } else object$CIobject <- NULL
  # MSEs computation
  if (inherits(object$fit,"HLfit")) {
    object$RMSEs <- .RMSEwrapper(object,verbose=FALSE)
  } else {
    if (eval_RMSEs) {
      if  (verbose) .overcat("Computing RMSEs... (may be slow)\n",prevmsglength)
      object$RMSEs <- .RMSEwrapper.SLik_j(object,verbose=FALSE)
    } else {
      object$pars <- .par_wrapper(object,verbose=FALSE) ## quick patch b/c no RMSEs
    } 
  }
  object$par_RMSEs <- .par_RMSEwrapper(object,verbose=FALSE)
  if  (verbose) {
    if ( ! is.null(object$par_RMSEs)) {
      .overcat("*** Interval estimates and RMSEs ***\n",prevmsglength)
      print(object$par_RMSEs)
    } else {
      .overcat("*** Interval estimates ***\n",prevmsglength)
      print(object$pars)
    }
  }
  if (is.null(object$logLs$cumul_iter)) {
    object$logLs$cumul_iter <- 1L ## adds a column to the data frame
    attr(object$logLs,"n_first_iter") <- nrow(object$logLs)
  }
  invisible(object)
}

