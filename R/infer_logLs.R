## FR->FR I should add the smoothing params to the return value of functions below (& for tailp)
##    so that I can perform smoothing and prediction on the smoothing params
##    NAs there should be allowed so that I don't resmooth predictions... complications expectable

infer_logL_by_GLMM <- function(EDF,stat.obs,logLname,verbose) {
#if (attr(EDF,"par")>0.97) browser()
  locEDF <- EDF[,names(stat.obs),drop=FALSE] 
  histo <- multi_binning(locEDF,focal=stat.obs)
  if (is.null(histo)) {
    #if (verbose) cat("!")
    obscov <- cov(locEDF) ## FR->FR remplacable par mixture model ?
    obsmean <- colMeans(locEDF)
    logL <- try(-LogAbsDetWrap(obscov) -log(2*pi)*length(obsmean) - (stat.obs-obsmean)%*%solve(obscov,(stat.obs-obsmean)),silent=TRUE)
    if (class(logL)=="try-error") { ## presumably because solve failed on singular matrix
      diag(obscov) <- diag(obscov)+ 1e-12 ## quick patch
      logL <- -LogAbsDetWrap(obscov) -log(2*pi)*length(obsmean) - (stat.obs-obsmean)%*%solve(obscov,(stat.obs-obsmean))
    }
    logL <- logL/2
    isValid <- FALSE
  } else {
    #if (verbose) cat(".")
    trend <- rep(1,nrow(histo)) 
    histo$trend <- trend
    fit <- smoothEDF_s(histo, pars=c())
    pred <- predict(fit,c(stat.obs,1,1)) ## 1 for binFactor and 1 for trend
    logL <- log(pred)
    isValid <- TRUE
  }
  names(logL) <- logLname
  unlist(c(attr(EDF,"par"),logL,isValid=isValid)) 
}


infer_logL_by_Hlscv.diag <- function(EDF,stat.obs,logLname,verbose) {
  locEDF <- EDF[,names(stat.obs),drop=FALSE] 
  chull <- resetCHull(locEDF[,names(stat.obs),drop=FALSE],formats=c("vertices","constraints")) 
  if ( ! isPointInCHull(stat.obs,constraints=chull[c("a","b")])) {
    #if (verbose) cat("!")
    ## quick gaussian approximation
    obscov <- cov(locEDF) ## FR->FR remplacable par mixture model ?
    obsmean <- colMeans(locEDF)
    logL <- (-determinant(obscov)$modulus -log(2*pi)*length(obsmean) 
             - (stat.obs-obsmean)%*%solve(obscov,(stat.obs-obsmean)) )
    logL <- logL/2
    isValid <- FALSE
  } else {
    #if (verbose) cat(".")
    if (requireNamespace("ks",quietly=TRUE)) {
      Hmat <- ks::Hlscv.diag(locEDF) ## bandwidth estimation
      fit <- ks::kde(locEDF,H=Hmat) ## smoothing
    } else {
      stop("ks package not available.")
    }
    logL <- log(predict(fit,x=stat.obs))    
    isValid <- TRUE
  }
  names(logL) <- logLname
  unlist(c(attr(EDF,"par"),logL,isValid=isValid))   
}

infer_logL_by_mclust <- function(EDF,stat.obs,logLname,verbose) {
  stats <- names(stat.obs)
  locEDF <- EDF[,stats,drop=FALSE] 
  if (requireNamespace("mclust",quietly=TRUE)) {
    fit <- mclust::densityMclust(locEDF)
  } else {
    stop("mclust package not available.")
  }
  pred <- predict(fit,t(c(stat.obs))) ## pas de binFactor!
  logL <- log(pred)
  if (invalid <- is.infinite(logL)) {
    #if (verbose) cat("!")
    logL <- -.Machine$double.xmax
  } #else if (verbose) cat(".")
  names(logL) <- logLname
  unlist(c(attr(EDF,"par"),logL,isValid= ! invalid)) 
}

infer_logL_by_Rmixmod <- function(EDF,stat.obs,logLname,verbose) {
  stats <- names(stat.obs)
  locEDF <- EDF[,stats,drop=FALSE] 
  fit <- densityMixmod(locEDF,stat.obs=stat.obs) ## Infusion::densityMixmod
  if (length(fit@nbCluster)==0L) { ## likely degenerate distribution
    checkfix <- apply(locEDF,2,var)==0
    ucheckfix <- apply(locEDF[,checkfix,drop=FALSE],2,unique)
    ## check that constant simulated statistics match stat.obs
    if (any(ucheckfix != stat.obs[checkfix])) {
      logL <- -.Machine$double.xmax
    } else {
      ## the likelihood can only be assessed by a continuity argument over parameter space
      invalid <- TRUE 
      logL <- NA
    }
  } else {
    logL <- predict(fit,tcstat.obs=t(c(stat.obs)),log=TRUE) ## pas de binFactor!
    if (invalid <- is.infinite(logL)) {
      #if (verbose) cat("!")
      logL <- -.Machine$double.xmax
    } #else if (verbose) cat(".")
  }
  names(logL) <- logLname
  unlist(c(attr(EDF,"par"),logL,isValid= ! invalid)) 
}

infer_logLs <- function(object,
                        stat.obs,
                        logLname=Infusion.getOption("logLname"), 
                        verbose=list(most=interactive(), ## must be explicitly set to FALSE in knitr examples
                                     final=FALSE),  
                        method="infer_logL_by_Rmixmod",
#                        fittedPars=NULL,
                        ... ## required because if generic includes them...
) {
  if ( ! is.null( cn <- colnames(stat.obs))) {
    message("Note: 'stat.obs' should be a numeric vector, not a matrix or data.frame.")
    names(stat.obs) <- cn ## minimal patch so that names() can be used, not colnames()
  }
  if (!is.list(verbose)) verbose <-as.list(verbose)
  if (is.null(names(verbose))) names(verbose) <- c("most","final")[seq_len(length(verbose))]
  if (is.null(verbose$most)) verbose$most <- interactive()
  if (is.null(verbose$final)) verbose$final <- FALSE
  allPars <- names(attr(object[[1]],"par"))
  object <- check_logLs_input(object)
  prevmsglength <- 0L
  nInvalid <- 0L
  lit <- 0L
  Sobs.densities <- lapply(object, function(element) {
    par_logL_indic <- do.call(method, list(EDF=element,stat.obs=stat.obs,logLname=logLname,verbose=verbose))
    # par_logL_indic is vector of pars + logL + isValid
    if (! par_logL_indic["isValid"]) nInvalid <<- nInvalid+1L
    lit <<- lit+1L
    if (verbose$most) {
      msg <- paste("Already ", lit, " distributions smoothed",sep="")
      if (nInvalid>0L) msg <- paste(msg," (",nInvalid," tagged as 'outlier'(s))",sep="")
      prevmsglength <<- overcat(msg, prevmsglength)
    }
    return(par_logL_indic)
  }) 
  if (verbose$final) {
    msg <- paste(lit, " distributions smoothed",sep="")
    if (nInvalid>0L) msg <- paste(msg," (",nInvalid," tagged as 'outlier'(s))",sep="")
    msg <- paste(msg,"            ")
    prevmsglength <<- overcat(msg, prevmsglength)
  }
  Sobs.densities <- do.call(rbind,Sobs.densities)  
  Sobs.densities <- data.frame(Sobs.densities[,c(allPars,logLname,"isValid"),drop=FALSE])
  if (nrow(Sobs.densities)==0L) {
    stop("No valid density estimates found")
  }
  attr(Sobs.densities,"EDFstat") <- method ## 
  attr(Sobs.densities,"stat.obs") <- stat.obs
  attr(Sobs.densities,"Simulate") <- attr(object,"Simulate")
  attr(Sobs.densities,"callfn") <- "infer_logLs"
  attr(Sobs.densities,"projectors") <- attr(object,"projectors")
  whichvar <- apply(Sobs.densities[,allPars,drop=FALSE],2L,function(v) length(unique(range(v)))>1)
  fixedPars <- Sobs.densities[1,!whichvar,drop=FALSE] ## has both names and values
  attr(Sobs.densities,"colTypes") <- list(allPars=allPars, ## keeps the order of the columns
                                          fittedPars=allPars[whichvar],
                                          fixedPars=fixedPars,
                                          logLname=logLname,
                                          statNames=names(stat.obs))
  attr(Sobs.densities,"LOWER") <- attr(object,"LOWER")
  attr(Sobs.densities,"UPPER") <- attr(object,"UPPER")
  attr(Sobs.densities,"Infusion.version") <- packageVersion("Infusion")
  
  class(Sobs.densities) <- c("logLs",class(Sobs.densities))
  return(Sobs.densities)
}

print.logLs <- function(x,...) {
  print.data.frame(x,...)
  addattrnames <- setdiff(names(attributes(x)),c("names","row.names","class"))
  cat(paste("Data frame with additional attributes:\n",paste(addattrnames,collapse=", "),"\n"))
}

summary.logLs <- function(object,...) {
  print(x=object,...)
}
