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
    lad <- determinant(obscov)$modulus[1]
    # pb general est cert eigenvalues peuvent -> +inf et d'autres -inf auquel cas logabsdet peut être innocuous mais pas estimaable précisément   
    if (is.nan(lad) || is.infinite(lad)){## because of determinant of nearly singular matrix
      zut <- abs(eigen(obscov,only.values = TRUE)$values) 
      zut[zut<1e-12] <- 1e-12
      lad <- sum(log(zut)) 
    }
    obsmean <- colMeans(locEDF)
    logL <- -lad -log(2*pi)*length(obsmean) - (stat.obs-obsmean)%*%solve(obscov,(stat.obs-obsmean))
    logL <- logL/2
    isValid <- FALSE
  } else {
    #if (verbose) cat(".")
    trend <- rep(1,nrow(histo)) 
    histo$trend <- trend
    fit <- .smoothEDF_s(histo, pars=c())
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
    Hmat <- .do_call_wrap("Hlscv.diag", list(x=locEDF), pack="ks")
    fit <- .do_call_wrap("kde", list(x=locEDF,H=Hmat), pack="ks")
    logL <- .do_call_wrap("predict", list(object=fit,x=stat.obs), pack="ks")
    isValid <- TRUE
  }
  names(logL) <- logLname
  unlist(c(attr(EDF,"par"),logL,isValid=isValid))   
}

infer_logL_by_mclust <- function(EDF,stat.obs,logLname,verbose) {
  stats <- names(stat.obs)
  locEDF <- EDF[,stats,drop=FALSE] 
  if ("package:mclust" %in% search()) { ## don't assume it was previous attached if so olny in a child process...
     ##requireNamespace + nclust:: not enough  given calls 
    fit <- ..densityMclust(locEDF,stat.obs=stat.obs) ## handling boundary effects, != mclust::densityMclust
  } else {
    stop("'mclust' should be loaded first.")
  }
  nbClu <- dim(fit$parameters$variance$sigma)[3L]
  solve_t_chol_sigma_list <- vector("list", nbClu)
  for (it in seq_len(nbClu)) solve_t_chol_sigma_list[[it]] <- solve(t(chol(fit$parameters$variance$sigma[,,it])))
  pred <- predict(fit,t(c(stat.obs)), solve_t_chol_sigma_list=solve_t_chol_sigma_list) ## pas de binFactor!
  logL <- log(pred)
  if (invalid <- is.infinite(logL)) {
    #if (verbose) cat("!")
    logL <- -.Machine$double.xmax
  } #else if (verbose) cat(".")
  names(logL) <- logLname
  unlist(c(attr(EDF,"par"),logL,isValid= ! invalid)) 
}

# Used for the primitive workflow
infer_logL_by_Rmixmod <- function(EDF,stat.obs,logLname,verbose) {
  stats <- names(stat.obs)
  locEDF <- EDF[,stats,drop=FALSE] 
  XtX <- crossprod(locEDF)
  if (is.infinite(kappa(XtX))) {
    warning("The summary statistics are collinear. Clustering will likely fail. Remove some summary statistic(s).")
  }
  fit <- .densityMixmod(locEDF,stat.obs=stat.obs) ## Infusion::densityMixmod # using function's default seed =Infusion.getOption("mixmodSeed")
  if (length(fit@nbCluster)==0L) { ## likely degenerate distribution
    checkfix <- sapply(locEDF,var)==0
    ucheckfix <- sapply(locEDF[,checkfix,drop=FALSE],unique)
    ## check that constant simulated statistics match stat.obs
    if (any(ucheckfix != stat.obs[checkfix])) {
      logL <- -.Machine$double.xmax
    } else {
      ## the likelihood can only be assessed by a continuity argument over parameter space
      invalid <- TRUE 
      logL <- NA
    }
  } else {
    solve_t_chol_sigma_list <- lapply(fit@parameters["variance"], .solve_t_cholfn)
    logL <- predict(fit,tcstat.obs=t(c(stat.obs)),solve_t_chol_sigma_list=solve_t_chol_sigma_list, 
                    clu_means=t(fit@parameters["mean",]),
                    logproportions = log(fit@parameters@proportions), log=TRUE) ## pas de binFactor!
    if (invalid <- is.infinite(logL)) {
      #if (verbose) cat("!")
      logL <- -.Machine$double.xmax
    } #else if (verbose) cat(".")
  }
  names(logL) <- logLname
  unlist(c(attr(EDF,"par"),logL,isValid= ! invalid)) 
}

## primitive workflow
## when the return value of infer_logLs() is changed the densv and densb data must be recomputed
infer_logLs <- function(object,
                        stat.obs,
                        logLname=Infusion.getOption("logLname"), 
                        verbose=list(most=interactive(), ## must be explicitly set to FALSE in knitr examples
                                     final=FALSE),  
                        method=Infusion.getOption("mixturing"),
                        nb_cores=NULL,
                        packages=NULL,
                        cluster_args=list(),
                        ... # currrently not used
) {
  if (is.data.frame(object)) {
    stop(paste0("'object' is a data.frame, not a list of matrices.\n ",
               "Did you mean to call infer_Slik_joint() rather than infer_logLs(),\n ",
               "or was Infusion's 'nRealizations' set to 1 globally?"))
  }
  if ( ! is.null( cn <- colnames(stat.obs))) {
    message("Note: 'stat.obs' should be a numeric vector, not a matrix or data.frame. Converting...")
    raw_data <- attr(stat.obs,"raw_data")
    stat.obs <- drop(stat.obs)
    attr(stat.obs,"raw_data") <- drop(raw_data)
    names(stat.obs) <- cn ## minimal patch so that names() can be used, not colnames()
  }
  if (!is.list(verbose)) verbose <- as.list(verbose)
  if (is.null(names(verbose))) names(verbose) <- c("most","final")[seq_len(length(verbose))]
  if (is.null(verbose$most)) verbose$most <- interactive()
  if (is.null(verbose$final)) verbose$final <- FALSE
  if (method=="mclust") method <- "infer_logL_by_mclust"
  if (method=="Rmixmod") method <- "infer_logL_by_Rmixmod"
  allPars <- names(attr(object[[1]],"par"))
  object <- .check_logLs_input(object)
  Sobs.densities <- vector("list", length(object)) 
  method_arglist <- list(stat.obs=stat.obs,logLname=logLname,verbose=verbose)
  #
  cluster_args <- .set_cluster_type(cluster_args, nb_cores) # PSOCK vs FORK
  cores_info <- .init_cores(cluster_args) ###### , ...)
  # make sure that a user-defined nondefault method is converted to a string: 
  if ( ! is.null(cores_info$cl) && ! is.null(nondefault <- match.call()$method)) method <- paste(nondefault) 
  if (method=="infer_logL_by_mclust") packages <- c(packages,"mclust")
  if (method=="infer_logL_by_Rmixmod") packages <- c(packages,"Rmixmod") ## did not appear necessary...
  Sobs.densities <- .eval_Sobs_densities(method, object, cores_info, stat.obs=stat.obs,packages=packages,logLname=logLname,verbose=verbose)
  .close_cores(cores_info)
  #
  Sobs.densities <- do.call(rbind,Sobs.densities)  
  Sobs.densities <- data.frame(Sobs.densities[,c(allPars,logLname,"isValid"),drop=FALSE])
  if (nrow(Sobs.densities)==0L) {
    stop("No valid density estimates found")
  }
  attr(Sobs.densities,"call") <- match.call() ## to get the extra arguments such as packages for refine() -> infer_logLs()
  attr(Sobs.densities,"EDFstat") <- method ## 
  attr(Sobs.densities,"stat.obs") <- stat.obs
  attr(Sobs.densities,"Simulate") <- attr(object,"Simulate")
  attr(Sobs.densities,"packages") <- list(add_simulation=attr(object,"packages"), logL_method=packages)
  attr(Sobs.densities,"env") <- attr(object,"env")
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
