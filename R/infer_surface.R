"infer_surface" <- function(object, ...) UseMethod("infer_surface") ## makes it easy to develop new inference methods

.check_surface_input <- function(object) {
  colTypes <- attr(object,"colTypes")
  fittedPars <- colTypes$fittedPars
  if (length(fittedPars)==0L) stop("No fitted parameters found.")
  checkX <- as.matrix(object[,fittedPars,drop=FALSE])
  checkX <- crossprod(checkX)
  checkX <- try(solve(checkX),silent=TRUE)
  if (inherits(checkX,"try-error")) {
    message(paste("Parameter names are",paste(fittedPars,collapse=", ")))
    stop("Some parameter(s) seem to be be a linear combination of other(s). Check input.")
  }
  return(TRUE)
} 

# object is a data frame with attributes; no object$colTypes or $stat.obs... nor $fit...
infer_surface.logLs <- function(object, method="REML",verbose=interactive(),allFix=NULL,...) {
  rownames(object) <- make.names(rownames(object),unique = TRUE)
  ## late addition 02/2016: 
  object <- object[ ! is.na(object[,"logL"]),,drop=FALSE]
  ## logL = NA then represents only singularities in distributions of statistics
  .check_surface_input(object)
  colTypes <- attr(object,"colTypes")
  fittedPars <- colTypes$fittedPars
  np <- length(fittedPars)
  logLname <- colTypes$logLname
  heuristic <- floor(50**((np/3)**(1/3))) ## 15  30  50  74 103 138 179 226 282 345
  if ( (nr <- nrow(object)) < heuristic) {
    message(paste("Only",nr,
                  "likelihood points available. It is advised to compute at least",
                  heuristic,"points."))
    heuristic <- min(nr,heuristic)
    quiteLow1 <- Inf ## to make sure that the second cluster will be eliminated (maybe not test quiteLow1 in any case? )
  } else quiteLow1 <- sort(object[,logLname],decreasing=TRUE)[heuristic]
  #quiteLow2 <- max(object[,logLname])-qchisq(1-1e-6,length(fittedPars))/2 ## ~ -24/2 -> -47/2 for 1 to 10 df
  quiteLow2 <- max(object[,logLname])-qchisq(1-0.1/Infusion.getOption("nRealizations"),length(fittedPars))/2 
  #quiteLow2 <- max(object[,logLname])-qchisq(1-0.001/Infusion.getOption("nRealizations"),length(fittedPars))/2 
  nrep <- nr -nrow(unique(object[,fittedPars,drop=FALSE]))
  if (nrep==0L) {
    stop("One *must* include replicate simulations in the 'object'.")
  }
  # ## try(,silent=TRUE) would be useless because clu@error will fail anyway. No error handling here:
  # clu <- suppressWarnings(mixmodCluster(object[,logLname],nbCluster=2,strategy = mixmodStrategy(seed=123))) 
  # if (clu@error) { ## Rmixmod has trouble fitting two clusters when single deviant point (zero-variance cluster)
  #   suspectPts <- which.min(object[,logLname])
  #   suspectPts <- (seq(nrow(object))==suspectPts) ## not elegant, but consistent with alternative code
  # } else {
  #   ## problem in identifying wrong points: when CIpoints accumulate, they form a low-variance cluster
  #   ## and all lower and higher points form another => we risk losing the higer points !
  #   ## => find the cluster of the lowest point
  #   suspectClu <- clu@results[[1L]]@partition[which.min(object[,logLname])]
  #   cluMeans <- clu@results[[1L]]@parameters@mean
  #   suspectPts <- (clu@results[[1L]]@partition == suspectClu & object[,logLname] < min(quiteLow1,quiteLow2)) ## may be all FALSE
  # }
  purgedlogLs <- object ## suspect points will be removed from the smoothing input but not from the return object's $logLs
  if (identical(Infusion.getOption("ESSAI"),TRUE)) {
    purgedlogLs[,logLname] <- pmax(purgedlogLs[,logLname],quiteLow2) - log(1-pmin(purgedlogLs[,logLname]-quiteLow2,0))
  } else {
    suspectPts <- try(.suspectPts_by_Rmixmod(logls=object[,logLname], threshold = min(quiteLow1,quiteLow2)),silent=TRUE)
    if (inherits(suspectPts,"try-error")) { ## if Rmixmod::mixmodCluster not accessible...
      suspectPts <- .suspectPts_by_mclust(logls=object[,logLname], threshold = min(quiteLow1,quiteLow2))
    }
    if (any(suspectPts)) {
      minGood <- min(purgedlogLs[ ! suspectPts,logLname])
      quiteLow <- minGood - min(1e-4,(minGood-max(purgedlogLs[suspectPts,logLname]))/1000) 
      purgedlogLs[suspectPts,logLname] <- NA ## tags them for matching and removal
      purgedlogLs <- .remove.pairswithNas(purgedlogLs) ## *removes pairs* where any member is NA   
      purgedlogLs <- purgedlogLs[ ! is.na(purgedlogLs[,logLname]),,drop=FALSE] ## removes *remaining NA's* not from pairs
    }  
  }
  #
  form <- as.formula(paste(logLname,"~ 1 + Matern(1|",paste(fittedPars,collapse=" + "),")"))
  lower <- apply(object,2,min)[fittedPars] ## this should be that of the full object as in the return $lower 
  upper <- apply(object,2,max)[fittedPars] ## idem
  if (is.null(allFix)) { ## estimation
    stat.obs <- attr(object,"stat.obs")
    if (verbose) cat(paste("\nUsing",method,"to infer the S-likelihood surface...\n"))
    if (method=="GCV") { ## the GCV rho estimates tend to be too low ?
      ## the effect of this GCV bias is over-smoothing, with phi not adjusted (since it comes from pure RMSE) 
      ## => overconfidence in precision of prediction; it may be better to have under confidence
      unLiklogLs <- prepareData(data=purgedlogLs,ParameterNames=fittedPars,
                                respName=logLname)
      gcvres <- calcGCV(unLiklogLs)
      ranfix <- list(rho=1/gcvres$CovFnParam[fittedPars],
                     nu=gcvres$CovFnParam[["smoothness"]])
      ## we need etaFix given ranFix and the same purgedlogLs
      thisfit <- HLCor(form,data=purgedlogLs,ranPars=ranfix,HLmethod="REML")
      ranfix <- c(ranfix,list(lambda=thisfit$lambda,phi=thisfit$phi))  
      etafix <- list(beta=fixef(thisfit))
    } else { ## handles all HLfit methods
      init.corrHLfit <- list(rho=1/(upper-lower))
      if (identical(Infusion.getOption("ESSAI"),TRUE)) {
        if (eval(Infusion.getOption("fitmeCondition"))) {
          thisfit <- fitme(form,data=purgedlogLs, fixed=list(nu=4),method=method, init=init.corrHLfit) 
        } else thisfit <- corrHLfit(form,data=purgedlogLs, ranFix=list(nu=4),
                                    HLmethod=method, init.corrHLfit=init.corrHLfit)
        # le code alternatif devrait marcher mais est complique et fragile
      } else {
        ## uses prior weights to reduce influence of extreme points but the prior weight computation is fragile.
        purgedlogLs$uli <- as.factor(.ULI(purgedlogLs[,fittedPars])) 
        init.phi.form <- as.formula(paste(logLname,"~ uli"))
        ## FR->FR palliatif a probleme de spaMM init values
        resglm <- glm(init.phi.form,data=purgedlogLs)
        init.phi <- as.numeric(deviance(resglm)/resglm$df.residual)
        init.phi <- min(9.99e7,max(1.01e-6,init.phi)) ## avoids extreme prior weights, but is ugly (F I X M E? here or in spaMM?)
        dlogL <- max(purgedlogLs$logL)-purgedlogLs$logL
        maxdlogL <- max(dlogL)
        priorwName <- .makenewname("priorw",names(purgedlogLs)) ## name not already in the data
        purgedlogLs[[priorwName]] <- (init.phi/sqrt(maxdlogL))^(dlogL/maxdlogL) ## notes 18/07/2016
        # then build a string without "priorwName" in it (but the value of priorwName) (2) convert to expression by parse() (3) eval
        if (eval(Infusion.getOption("fitmeCondition"))) {
          thisfit <- eval(parse(text=paste("fitme(form,data=purgedlogLs, fixed=list(nu=4),", 
                                           "method=method, init=init.corrHLfit,prior.weights=",priorwName,")"))) 
        } else thisfit <- eval(parse(text=paste("corrHLfit(form,data=purgedlogLs, ranFix=list(nu=4),init.HLfit=list(phi=init.phi),", 
                                                "HLmethod=method, init.corrHLfit=init.corrHLfit,prior.weights=",priorwName,")")))
      }
      RMSE <- sqrt(thisfit$phi)
      if (thisfit$spaMM.version<"2.4.26") {
        corrPars1 <- thisfit$corrPars[["1"]]
      } else corrPars1 <- get_ranPars(thisfit,which="corrPars")[["1"]]
      ranfix <- c(corrPars1,list(lambda=thisfit$lambda,phi=thisfit$phi))
      ## it is worth fixing the fixed effects if all other params are fixed. Otherwise 
      ## extreme low response values in the 'object' affect the predictionCoeffs and create artificial "valleys" in the predictions
      ## (it might be better to refit lambda, except that we don't want extreme lambda...)  
      etafix <- list(beta=fixef(thisfit)) 
      if (FALSE) { ## relation between the two rho estimates
        unLiklogLs <- prepareData(data=purgedlogLs,ParameterNames=fittedPars,
                                  respName=logLname)
        blackbox.options(maxSmoothness=4) # !
        gcvres <- calcGCV(unLiklogLs)
        print(c(ranfix$rho,1/gcvres$CovFnParam[fittedPars]))
      }
    } 
  } else if (method=="newdata") { ## keep allFix, produces new predictor from new points
    boundedlogLs <- object   
    if (any(suspectPts)) boundedlogLs[suspectPts,logLname] <- quiteLow
    etafix <- allFix["beta"]
    if (is.null(etafix)) stop("'etafix' is NULL for method = newdata in infer_surface.logLs")
    ranfix <- allFix
    ranfix$beta <- NULL
    # With etaFix vs without it, p_bv is computed differently, except if REMLformula is used
    thisfit <- corrHLfit(form,data=boundedlogLs,ranFix=ranfix,etaFix=etafix,REMLformula=form)
#    if (verbose) cat(paste("for newdata:",thisfit$APHLs$p_bv,"\n"))    
  } else { ## 'allFix' is TRUE and 'method' is not "newdata"
    etafix <- allFix["beta"]
    ranfix <- allFix
    ranfix$beta <- NULL
    # With etaFix vs without it, p_bv is computed differently, except if REMLformula is used
    thisfit <- corrHLfit(form,data=purgedlogLs,ranFix=ranfix,etaFix=etafix,REMLformula=form) 
#    if (verbose) cat(paste("for testing:",thisfit$APHLs$p_bv,"\n"))    
  }
  # 
  .Object <- list()
  .Object$rho <- thisfit$ranFix$rho # predict.Density uses this
  .Object$colTypes <- colTypes
  .Object$corr.args <- thisfit$ranFix[which(names(thisfit$ranFix) %in% c("nu","Nugget"))]
  ###  .Object$offset.formula <- thisoffset
  ## FR->FR should ultimately use convex hulls ?
  .Object$lower <- lower
  .Object$upper <- upper
  if (any(is.na(c(lower,upper)))) stop("NA(s) in c(lower,upper))")
  .Object$LOWER <- attr(object,"LOWER")
  .Object$UPPER <- attr(object,"UPPER")
  # thisfit$predictionCoeffs <- predictionCoeffs(thisfit)
  .Object$logLs <- object ## (fittedPars,logL) and attributes (stat.obs) Simulate packages env!
  obspred <- predict(thisfit,variances=list(linPred=TRUE,dispVar=TRUE),binding=logLname)
  .Object$obspred <- obspred
  # Qmax
  obsSE <- attr(obspred,"predVar")
  obsSE[obsSE<0] <- 0 ## anticipating numerical problems
  .Object$Qmax <- max(obspred[,attr(obspred,"fittedName")]+1.96 * sqrt(obsSE)) ## best improvement function for already computed points 
  #
  .Object$fit <- thisfit   
  .Object$projectors <- attr(object,"projectors")
  .Object$EDFstat <- attr(object,"EDFstat")## may be NULL
  .Object$`Infusion.version` <- packageVersion("Infusion")
  #.Object$`ori.stats.names` <- attr(object,"ori.stats.names")
  class(.Object) <- c("SLik",class(.Object))
  return(.Object)   
} 

## version that overcome new problems in Rmixmod 2.1.0: we need one clustering to succeed.
# 1 cluster should always work, so nbCluster=1:2 even if only2 is interesting 
.suspectPts_by_Rmixmod <- function(logls,threshold) {
  clu <- suppressWarnings(.do_call_wrap("mixmodCluster",list(data=logls,nbCluster=1:2,seed=123) ))
  if (clu@results[[2L]]@error=="No error") {
    ## problem in identifying wrong points: when CIpoints accumulate, they form a low-variance cluster
    ## and all lower and higher points form another => we risk losing the higer points !
    ## => find the cluster of the lowest point
    suspectClu <- clu@results[[2L]]@partition[which.min(logls)]
    suspectPts <- (clu@results[[2L]]@partition == suspectClu & logls < threshold) ## may be all FALSE
  } else {
    suspectPts <- which.min(logls)
    suspectPts <- (seq(length(logls))==suspectPts) ## not elegant, but consistent with alternative code
  }
}

.suspectPts_by_mclust <- function(logls,threshold) {
  if ("package:mclust" %in% search()) { 
    clu <- .do_call_wrap("Mclust",
                         list(data=logls,G=2,verbose=FALSE),
                         pack="mclust")
  } else  stop("'mclust' should be loaded first.") ## occurs if only loaded in a child process...
  ## problem in identifying wrong points: when CIpoints accumulate, they form a low-variance cluster
  ## and all lower and higher points form another => we risk losing the higer points !
  ## => find the cluster of the lowest point
  suspectClu <- clu$classification[which.min(logls)]
  suspectPts <- (clu$classification == suspectClu & logls < threshold) ## may be all FALSE
}

