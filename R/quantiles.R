.infer_LR_pvalue_from_EDF  <- function(EDF,stat.obs,tailNames,
                             refDensity, ## a reference density used for ordering of all EDF samples; typically the density( ; ML)
                             verbose) {
  locEDF <- EDF[,names(stat.obs),drop=FALSE] 
  histo <- multi_binning(locEDF,focal=stat.obs)
  if (is.null(histo)) {
    if (verbose) cat("!")
    tailp <- 0 
    isValid <- FALSE
  } else {
    if (verbose) cat(".")
    trend <- rep(1,nrow(histo)) 
    histo$trend <- trend
    fit <- .smoothEDF_s(histo, pars=c())
    L0EDF <- as.numeric(predict(fit,cbind(EDF,binFactor=1,trend=1))) 
    L1EDF <- as.numeric(predict(refDensity,cbind(EDF,binFactor=1,trend=1))) 
    L0obs <- as.numeric(predict(fit,c(stat.obs,1,1))) ## 1 for binFactor and 1 for trend
    L1obs <- as.numeric(predict(refDensity,c(stat.obs,1,1))) ## 1 for binFactor and 1 for trend
    # EDF being simulated ~H0, ref ~H1 : p-value= 
    tailp <- length(which((L0EDF/L1EDF<(L0obs/L1obs)))) 
    isValid <- TRUE
  }
  tailp <- c(tailp,nrow(EDF)-tailp)
  names(tailp) <- tailNames
  unlist(c(attr(EDF,"par"),tailp,isValid=isValid)) 
}
# wrongordering.R.txt has a more "Fisherian" ordering function

# fixme not currently used
.infer_refDensity <- function(object,nRealizations= ## revise
                               10*Infusion.getOption("nRealizations")) { ## from SLik object
  ori_nReal <- Infusion.getOption("nRealizations")
  refSimul <- add_simulation(,Simulate=get_from(object,"Simulate"),
                             nRealizations = nRealizations,
                             parsTable=t(c(object$MSL$MSLE,object$colTypes$fixedPars)))
  Infusion.options(nRealizations=ori_nReal)
  stat.obs <- get_from(object,"stat.obs")
  EDF <- (refSimul[[1]])[,names(stat.obs),drop=FALSE]
  histo <- multi_binning(EDF,focal=stat.obs)
  histo$trend <- 1
  fit <- .smoothEDF_s(histo, pars=c())
  return(fit)
}


infer_tailp <- function(object, ## object is a list of EDFs
                            refDensity,
                            stat.obs,
                            tailNames=Infusion.getOption("tailNames"), 
                            verbose=interactive(), method=NULL, cluster_args, 
                            ... ## required because if generic includes them...
) {
  
  if (is.null(method)) {
    locfn <- ".infer_LR_pvalue_from_EDF"
  } else locfn <- method ## leaves method unchanged, see attributes of return value below   
  allPars <- names(attr(object[[1]],"par"))
  Sobs.tailp <- lapply(object, locfn, refDensity=refDensity, stat.obs=stat.obs,tailNames=tailNames,verbose=verbose) 
  Sobs.tailp <- do.call(rbind,Sobs.tailp)  
  Sobs.tailp <- data.frame(Sobs.tailp[,c(allPars,tailNames,"isValid"),drop=FALSE])
  ## Attributs encombrants pour une data.frame. Puisqu'on a défini une classe, on pourrait la structurer autrement... 
  ##           le pb c'est celui de traitements homogènes avec SLik 
  attr(Sobs.tailp,"call") <- match.call()
  attr(Sobs.tailp,"EDFstat") <- method ## maybe not useful see in "call" : fixme
  attr(Sobs.tailp,"stat.obs") <- stat.obs ## maybe not useful see in "call" : fixme
  attr(Sobs.tailp,"refDensity") <- refDensity ## maybe not useful see in "call" : fixme
  attr(Sobs.tailp,"Simulate") <- attr(object,"Simulate")
  attr(Sobs.tailp,"callfn") <- "infer_tailp" ## maybe not useful see in "call" : fixme
  attr(Sobs.tailp,"projectors") <- attr(object,"projectors")
  whichvar <- apply(Sobs.tailp[,allPars,drop=FALSE],2,function(v) length(unique(range(v)))>1)
  fittedPars <- allPars[whichvar]
  fixedPars <- Sobs.tailp[1, ! whichvar,drop=FALSE] ## has both names and values
  attr(Sobs.tailp,"colTypes") <- list(allPars=allPars, ## keeps the order ot he columns
                                          fittedPars=fittedPars,
                                          fixedPars=fixedPars,
                                          tailNames=tailNames,
                                          statNames=names(stat.obs)) 
  class(Sobs.tailp) <- c(class(Sobs.tailp),"tailp")
  return(Sobs.tailp)
}




confint.SLikp <- function(object, parm, ## parm is the parameter which CI is sought 
                           level=0.95, verbose=interactive(),fixed=NULL,which=c(TRUE,TRUE),...) {
  # conversion from LRT n df to 1 df to have a 1D confidence interval
  level <- pchisq( qchisq(level, df=length(object$MSL$MSLE)), df=1) ## e.g. 0.9856 pour 2df
  .confintAll(object=object, parm=parm, ## parm is the parameter which CI is sought 
             givenmax = 1,
             ecdf_2lr=NULL,
             level= - level, verbose=verbose,fixed=fixed,which=which,...)
}

# moved from spaMM

.noNonSpatialbarsMM <- function (term) {
  if (!("|" %in% all.names(term))) 
    return(term)
  if (is.call(term) && term[[1]] == as.name("|")) 
    return(NULL) ## removes (|) but not Matern(|)
  if (length(term) == 2L) {
    term1 <- as.character(term[[1]])
    if (term1 %in% c("adjacency","Matern","AR1","corrMatrix","ar1")) {
      return(term) 
    } else return(.noNonSpatialbarsMM(term[[2]])) 
  }
  nb2 <- .noNonSpatialbarsMM(term[[2]])
  nb3 <- .noNonSpatialbarsMM(term[[3]])
  if (is.null(nb2)) 
    return(nb3)
  if (is.null(nb3)) 
    return(nb2)
  term[[2]] <- nb2
  term[[3]] <- nb3
  term
}


infer_surface.tailp <- function(object, 
                                    method="PQL", ## default changed 23/03/2015 
                                    verbose=interactive(),
                                    allFix=NULL,
                                    ...) {
  .check_surface_input(object)
  colTypes <- attr(object,"colTypes")
  fittedPars <- colTypes$fittedPars
  EDFestLevelName <- Infusion.getOption("EDFestLevelName")
  .Object <- list()
  .Object$tailp <- object ## before any modification (bc addition of the <EDFestLevelName> would make rbinding with new tailp complicated)
  .Object$refDensity <- attr(object,"refDensity")
  attr(object,"refDensity") <- NULL ## so that it is not stored in the fitme outputs...!  
  object[,EDFestLevelName] <-  - seq(nrow(object)) ## '-' makes it easy to define new levels in 'predict'
  ## : (1|<EDFestLevelName>) needed to represent the extra variability of the density estimation procedure (much better)
  purgedlogLs <- object
  tailNames <- colTypes$tailNames
  purgedlogLs <- purgedlogLs[ purgedlogLs[,tailNames[1]]>0,,drop=FALSE] 
  form <- paste(tailNames,collapse=",")
  form <- as.formula(paste("cbind(",form,") ~ 1 + Matern(1|",paste(fittedPars,collapse=" + "),") + (1|",EDFestLevelName,")"))
  lower <- sapply(object,min)[fittedPars] ## this should be that of the full object as in the return $lower 
  upper <- sapply(object,max)[fittedPars] ## idem
  if (is.null(allFix)) {
    init <- list(rho=1/(upper-lower))
    if  (method=="GCV") {
      stop("GCV method not meaningful for infer_surface.tailp. Use PQL instead.")
      method <- "PQL" ## possible fallback 
    }
    if (verbose) cat(paste("\nusing",method,"to infer the S-tail p surface...\n"))
    thisfit <- fitme(form,data=purgedlogLs,
                         family=binomial(),
                         fixed=list(nu=4),
                         HLmethod=method,
                         init=init)
    corrPars1 <- get_ranPars(thisfit,which="corrPars")[["1"]]
    allFix <- c(corrPars1,list(lambda=thisfit$lambda))
  } else {
  }
  # 
  thisfit <- fitme(form,data=object,family=binomial(),HLmethod=method,fixed=allFix) ## full data smoothed with ranFix from cleaned data
  .Object$rho <- thisfit$ranFix$rho # density computation uses this
  .Object$colTypes <- colTypes
  .Object$corr.args <- thisfit$ranFix[which(names(thisfit$ranFix) %in% c("nu","Nugget"))]
  .Object$lower <- lower
  .Object$upper <- upper
  # thisfit$predictionCoeffs <- predictionCoeffs(thisfit)
  obspred <- predict(thisfit,variances=list(predVar=TRUE),binding=tailNames[1]) 
  .Object$obspred <- obspred ## here is why binding is used
  # Qmax
  obsVar <- attr(obspred,"predVar")
  if (any(obsVar<0))  obsVar <- 0 ## anticipating numerical problems (ignore also spuriously large values) ## but then Qmax will be ignored...
  #dmudeta <- thisfit$family$mu.eta(obspred[,tailNames[1]]) 
  #obsVar <- obsvar * dmudeta^2 
  #.Object$Qmax <- max(thisfit$family$linkfun(obspred[,attr(obspred,"fittedName")])+1.96 * sqrt(obsVar)) ## best improvement function for computed points
  obspred <- obspred[,attr(obspred,"fittedName")] ## unbinds
  # EI on linear predictor scale vs prediction on response scale
  .Object$Qmax <- max(thisfit$family$linkfun(obspred)+1.96 * sqrt(obsVar)) ## best improvement function for computed points 
  #
  .Object$fit <- thisfit   
  .Object$predictReForm <- .noNonSpatialbarsMM(form)
  .Object$projectors <- attr(object,"projectors")
  .Object$EDFstat <- attr(object,"EDFstat")## may be NULL
  class(.Object) <- c("SLikp",class(.Object))
  return(.Object)   
} 

refine.SLikp <- function(object,method=NULL,...) {
  if (is.null(method)) method <- "PQL"
  refine.default(object,surfaceData=object$tailp,method=method,...)
}

print.SLikp <-function(x,...) {summary(x)} ## indeed summarizes the list... unless a summary.SLik definition is available

`summary.SLikp` <- function(object,...) {
  if ( !is.null(object$MSL) ) {
    cat("*** Summary MQ: ***\n")
    print(c(object$MSL$MSLE,"tailp"=object$MSL$maxlogL))
    #
    if( ! is.null(CIobject <- object$CIobject))  {
      if (is.null(wrn <- object$CIobject$warn)) {
        cat("*** Confidence intervals: ***\n")
        cis <- lapply(CIobject$CIs,function(lt) {lt$interval})
        print(cis)
        locst <-  "*** RMSEs of MaxSumm-tail p and of MSQ for CIs: ***\n"
      } else {
        message(wrn)
        locst <-  "*** RMSE of MaxSumm-tail p: ***\n"    
      }
    } else {
      locst <-  "*** RMSE of MaxSumm-tail p: ***\n"    
    }
    cat(locst)
    if (is.null(wrn <- object$RMSEs$warn)) {
      print(get_from(object,"RMSEs"))
    } else message(wrn)
  } else {
    cat("SLikp object created. Use MSL(.) to obtain point estimates and CIs.\n")
  }
  invisible(object)
}

calc.lrthreshold.SLikp <- function(object,dlr=NULL,verbose=interactive(),...) {
  if (is.null(dlr)) dlr <- Infusion.getOption("LRthreshold")
  probErr <- object$MSL$predVar *2.326348 # qnorm(0.99,0,1) ## $predVar differs in conception from the pred MSE used in migraine
  ## convert to tailp then logit tailp
  dtailp <- pchisq(-dlr*2,df=1) ## 0.999 si LRthreshold= - qchisq(0.999,df=1)/2
  dlogit <- object$fit$family$linkfun(1-dtailp) ## typically logit(0.001)
  if (verbose) {
    cat("Default logit(tailp) threshold and probable prediction error:\n")
    locstring <- paste(.prettysignif(-dlr)," and ",.prettysignif(probErr),"\n",sep="")
    cat(locstring)
  }
  ## expand beyond *predicted* dlr threshold  ## this should be disconnected from GV$hullExpandFactor
  return( object$fit$family$linkinv(dlogit*1.2 -probErr) ) ## returns a tail p threshold
}

predict.SLikp <- function (object, newdata=object$tailp[,object$colTypes$fittedPars,drop=FALSE],...) {
  if (is.vector(newdata)) {
    newdata <- c(newdata,1) ## adds an <EDFestLevelName> value different from any preexisting one (all <0)
  } else newdata[,Infusion.getOption("EDFestLevelName")] <- seq(nrow(newdata)) ## idem
  return(predict(object$fit,newdata=newdata,...))
}

plot.SLikp <-function(x,y,filled=FALSE,log.=TRUE,...) {
  object <- x
  fittedPars <- object$colTypes$fittedPars
  resp <- object$obspred[,attr(object$obspred,"fittedName")] 
  saferespmin <- min(resp[resp>0])/2 
  if (log.) {
    Ztransf <- function(Z) log(pmax(Z,saferespmin))
    mainst <- "log(summary tail p) surface"
  } else {
    Ztransf <- NULL
    mainst <- "Summary tail p surface"
  }
  if (length(fittedPars)==1L) {
    x <- object$obspred[,fittedPars]
    if (log.) resp <- Ztransf(resp) 
    plot(x,resp,main="log(summary tail p)",xlab=fittedPars,ylab="log(S-tail p)")
  } else if (length(fittedPars)==2L) {
    if (inherits(object$fit,"HLfit")) {
      if (filled) {
        filled.mapMM(object$fit,
                     Ztransf=Ztransf,
                     color.palette=.Inf_palette(variant="spaMM_shift3"),nlevels=50,
                     plot.title=quote(title(main=mainst,
                                            xlab=fittedPars[1],ylab=fittedPars[2])),
                     decorations=quote({if (!is.null(object$latestPoints)) points(object$tailp[object$latestPoints,fittedPars],pch=".",cex=2);
                                        if (!is.null(object$CIobject)) points(object$CIobject$bounds,pch="+");
                                        points(t(object$MSL$MSLE),pch="+");
                                        points(object$tailp[,fittedPars],cex=0.8);
                                        points(object$tailp[ ! object$tailp[,"isValid"],fittedPars],pch=20,cex=0.8)}),...)
        
      } else mapMM(object$fit,
                   Ztransf=Ztransf,
                   color.palette=.Inf_palette(variant="spaMM_shift3"),nlevels=50,
                   plot.title=quote(title(main=mainst,
                                          xlab=fittedPars[1],ylab=fittedPars[2])),
                   decorations=quote({if (!is.null(object$latestPoints)) points(object$tailp[object$latestPoints,fittedPars],pch=".",cex=2);
                                      if (!is.null(object$CIobject)) points(object$CIobject$bounds,pch="+");
                                      points(t(object$MSL$MSLE),pch="+");
                                      points(object$tailp[ ! object$tailp[,"isValid"],fittedPars],pch=20,cex=0.8)}),...)
    } 
  } else if (length(fittedPars)>2L) {message("No plot.SLikp output for length(fittedPars)>2")}
  invisible(object)
}
