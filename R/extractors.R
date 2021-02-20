"get_from" <- function(object, which, ...) UseMethod("get_from") 

`get_from.default` <- function(object, which, ...) {
  return(object[[which]])
}
  
get_from.SLik <- function(object, which, raw=FALSE, ...) {
  if (raw) return(object[[which]]) # per the doc.
  ###  all (! raw) cases
  ## Specific 'which' for old versions 
  if (object$Infusion.version > "1.4.13") {
    if (which=="MSL") {
      return(as.list(object[["MSL"]]))
    } else if (which=="par_RMSEs") {
      return(object[["par_RMSEs"]]$par_RMSEs)
    } else if (which=="RMSEs") {
      return(object[["RMSEs"]]$RMSEs)
    }
  }
  # for all versions
  if (which %in% c("obs","stat.obs"))  return(attr(object$logLs,"stat.obs")) 
  # All other cases:
  return(object[[which]])
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
    interval <- CIs[[parm]]$interval
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
