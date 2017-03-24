
calc_EI <- function(summInferObject,points,Qmax=NULL) {  ## for both SLik and SLikp
  trypred <- predict(summInferObject,newdata=points,variances=list(linPred=TRUE,dispVar=TRUE))
  tryVar <- attr(trypred,"predVar")
  if (any(tryVar<0))  { ## anticipating numerical problems (ignore also spuriously large tryVar values)
    return( summInferObject$fit$family$linkfun(as.numeric(trypred[,1])) )
  } else {
    trySE <- sqrt(tryVar)
    #dmudeta <- summInferObject$fit$family$mu.eta(as.numeric(trypred)) ## as.numeric also removes attributes, quite useful
    #tryVar <- tryVar * dmudeta^2  
    ## EI on linear predictor scale vs prediction on response scale
    tryQ <- summInferObject$fit$family$linkfun(as.numeric(trypred[,1])) + 1.96*trySE ## improvement function for candidate points
    if (is.null(Qmax)) Qmax <- max(tryQ)
    return( trySE*dnorm((Qmax-tryQ)/trySE)+(tryQ-Qmax)*pnorm((tryQ-Qmax)/trySE) )## 7.5 p. 121              
  }
}

"calc.lrthreshold" <- function(object, ...) UseMethod("calc.lrthreshold") ## makes it easy to develop new inference methods

# has an SLik method

calc.lrthreshold.default <- function(object,dlr=NULL,verbose=interactive(),...) {
  stop(pastefrom("No default method for calc.lrthreshold. Provide a method."))
}


"refine" <- function(object, ...) UseMethod("refine") ## makes it easy to develop new inference methods

# has SLikp and SLik method

## si Simulate est exterieure, il faut que l'utilisateur puisse decomposer la fn et sample_volume doit être public...
`refine.default` <- function(object, ## SLikp or SLik
                             surfaceData, ## object$logLs or object$tailp, with stat.obs attribute, etc
                             maxit=1,n=NULL,useEI = list(max=TRUE,profileCI=TRUE,rawCI=FALSE),newsimuls=NULL,
                             useCI=TRUE,level=0.95,verbose=list(most=interactive(),movie=FALSE),
                             precision = Infusion.getOption("precision"),
                             method, ## "GCV" and HLfit methods for Slik objects
                             ...) {
  if (!is.list(verbose)) verbose <- as.list(verbose)
  if (is.null(names(verbose))) names(verbose) <- c("most","movie")[seq_len(length(verbose))]
  if (is.null(verbose$most)) verbose$most <- interactive()
  if (is.null(verbose$movie)) verbose$movie <- FALSE
  it <- 0L
  previous_cumul_iter <- max(object$logLs$cumul_iter)
  RMSEs <- object$RMSEs
  if (is.null(RMSEs)) RMSEs <- 1e10
  stat.obs <- attr(surfaceData,"stat.obs")
  if ( ! is.null(newsimuls) && maxit>1) stop("'maxit'>1 is incompatible with user-provided 'newsimuls'") 
  if  ( ! any(is.na(RMSEs)) && all(RMSEs<precision)) {
    cat("Target precision appears to be already reached.\n")
    # return(object) 
    ## nevertheless continue for one iteration 
  }
  EIsampleFactor <- 30
  while( it==0L ##always perform one iteration on request  
         || (it <maxit && (any(RMSEs>precision) || any(is.na(RMSEs)))) ) {
    if (maxit>1L && verbose$most) cat(paste("iter = ",it+1L,":\n",sep=""))
    fittedPars <- object$colTypes$fittedPars
    logLname <- object$colTypes$logLname
    if ( is.null(newsimuls)) {
      if (inherits(object,"SLik_j")) {
        locblob <- .rparam_from_SLik_j(object,fittedPars=fittedPars,level=level) 
        trypoints <- locblob$trypoints
        freq_good <- locblob$freq_good
      } else {
        if (verbose$most) cat("\nDesigning new parameter points (may be slow)...\n")
        ## replace pairs with low predicted lik by pairs with high predicted lik
        pred <- predict(object$fit,binding=logLname) 
        ## : corrected 11/07/2016: pred was predict(object,.) using object$logLs, not object$fit, thereby not matching uli
        uli <- ULI(object$fit$data[,fittedPars])
        table_uli <- table(uli)
        CIthreshold <- object$MSL$maxlogL - (qchisq(level,df=1)/2)
        singletsLvls <- names(table_uli)[table_uli==1L]
        singletsBool <- uli %in% singletsLvls
        single_pred <- pred[singletsBool,] ## (non-redundant) table of predicted values
        sort_single_pred <- single_pred[order(single_pred[,logLname],decreasing=TRUE),]
        doubletsBool <- ! uli %in% singletsLvls
        double_pred <- pred[doubletsBool,]
        unique_pred <- unique(double_pred)
        sort_unique_pred <- unique_pred[order(unique_pred[,logLname],decreasing=FALSE),]
        ncomp <- min(nrow(sort_unique_pred),nrow(sort_single_pred))
        if (ncomp>0L) {
          n_sub <- max(0,which(sort_single_pred[1:ncomp,logLname] > sort_unique_pred[1:ncomp,logLname]))
        } else n_sub <- 0L ## number than can be substituted
        if (n_sub>0L) {
          n_top <- length(which(sort_single_pred[,logLname]> CIthreshold))
          if (n_top< n_sub) {
            goodrows <- seq(n_sub) ## takes all highest
          } else goodrows <- sample(n_top,n_sub) ## random sampling of top  
          goodrows <- rownames(sort_single_pred[goodrows,,drop=FALSE])
          newpairs <- surfaceData[goodrows,object$colTypes$allPars,drop=FALSE]
          ## REMOVE ONE REPLICATE OF EACH "POOR" PAIR IN object$logLs
          surfaceData <- surfaceData[ ! rownames(surfaceData) %in% rownames(sort_unique_pred)[seq(n_sub)],]
        } else newpairs <- NULL
        # and the really slow part: 
        trypoints <- do.call(Infusion.getOption("rparamfn"),
                             list(object=object,n=n,useEI=useEI,tryn=EIsampleFactor*n,useCI=useCI,level=level,verbose=FALSE))
        trypoints <- rbind(trypoints,newpairs)
      }
      Simulate <- attr(surfaceData,"Simulate")
      if (is.null(Simulate)) return(trypoints)
      if (inherits(object,"SLik_j")) {
        newsimuls <- add_reftable(Simulate=Simulate,par.grid=trypoints,verbose=verbose$most)     
        newsimuls$cumul_iter <- previous_cumul_iter + it +1L
      } else {
        newsimuls <- add_simulation(Simulate=Simulate,par.grid=trypoints,verbose=verbose$most)    
      }
    }
    projectors <- object$projectors
    if ( ! is.null(projectors)) newsimuls <- project(newsimuls,projectors=eval(projectors)) 
    if (inherits(object,"SLik_j")) {
      if ( ! (is.matrix(newsimuls) || is.data.frame(newsimuls)) ) {
        stop("'newsimuls' must be a matrix or data.frame for refine.Slik_j() method.")
      }
      checkNA <- apply(newsimuls,1L,anyNA)
      newsimuls <- newsimuls[!checkNA,,drop=FALSE] ## FR->FR quick patch
      jointEDF <- structure(rbind(object$logLs,newsimuls), allPars=object$colTypes$allPars)
      object <- infer_SLik_joint(data = jointEDF,
                            stat.obs = attr(object$logLs,"stat.obs"),
                            Simulate = attr(object$logLs,"Simulate"),
                            verbose=verbose$most)
      object$latestPoints <- nrow(jointEDF)+1-seq_len(nrow(newsimuls)) ## for plots
    } else {
      if ( ! inherits(newsimuls,"list") ) {
        stop("'newsimuls' must be a _list_ of empirical distributions.")
      }
      arglist <- list(object=newsimuls,stat.obs=stat.obs,verbose=verbose$most,method=object$EDFstat) ## method may be NULL
      if (inherits(object,"SLikp")) arglist$refDensity <- object$refDensity 
      newlogLs <- do.call(attr(surfaceData,"callfn"),arglist)   ### call density estimation
      if (verbose$most) cat ("\n")
      successrate <- length(which(newlogLs$isValid>0))/nrow(newlogLs)
      EIsampleFactor <- EIsampleFactor * 0.98/successrate
      newlogLs$cumul_iter <- previous_cumul_iter + it +1L
      surfaceData <- rbind(surfaceData,newlogLs)
      itmethod <- method[min(length(method),it+1L)] ## may be overriden below when hat(nu) is low.
      ## tests whether resmoothing can yield substantial improvements:
      allFix <- c(object$fit$corrPars[c("nu","rho")],list(lambda=object$fit$lambda[1],phi=object$fit$phi[1],beta=fixef(object$fit)))
      previousRho <- allFix$rho ## full length vector of scale params
      if (is.null(previousRho)) {
        stop("is.null(previousRho) in 'refine.default'. Check code (trRho?).")
      }
      smoothingOK <- (allFix$nu>3.95) ## provisional assessment
      if (smoothingOK) { ## further assessments
        # Test whether the old fit predicts well the new points
        prednewfromold <- predict(object,newdata=newlogLs[,fittedPars,drop=FALSE],
                variances=list(linPred=TRUE,disp=TRUE))  
        msepred <- mean((prednewfromold[,1]-newlogLs[,logLname])^2)
        respvar <- mean(attr(prednewfromold,"predVar")) + object$fit$phi ## FR->FR does not account forprior weights
        relerr <- respvar/(1e-6+msepred)
        #print(paste("relerr:",relerr))
        smoothingOK <- (relerr>0.8 && relerr<1.25) ## practically always, given distrib of estimator 'relerr'  
        #
        currsurf <- infer_surface(surfaceData,method=itmethod,verbose=FALSE,allFix=allFix)
        #
        if (FALSE) {
          # assessment by perturbing rho
          currp_bv <- currsurf$fit$APHLs$p_bv
          testCorr <- allFix
          testCorr$rho <- previousRho*1.1
          newp_bv <- infer_surface(surfaceData,method=itmethod,verbose=FALSE,allFix=testCorr)$fit$APHLs$p_bv
          if ( smoothingOK <- (newp_bv < currp_bv) ) { ## no progress by increasing rho
            testCorr$rho <- previousRho/1.1
            newp_bv <- infer_surface(surfaceData,method=itmethod,verbose=FALSE,allFix=testCorr)$fit$APHLs$p_bv
            smoothingOK <- (newp_bv < currp_bv) ## no progress by decreasing rho
          }
        }
      } else if (itmethod=="GCV") itmethod <- "REML" ## (overrides GCV user option) enforces REML when GCV does not estimate nu correctly
      if (smoothingOK) {
        #print(paste("smoothingOK:",smoothingOK))
        # object <- infer_surface(surfaceData,method="newdata",verbose=FALSE,allFix=allFix) 
        ## newdata => source des biais évidents: pas assez de pointsenhaut de la surface ?
        object <- currsurf ## ie infer surface from purgedlogLs from new data,method=itmethod, with allFix=<previous pars>
      } else {
        object <- infer_surface(surfaceData,method=itmethod,verbose=verbose$most) ## new smoothing;could reuse corrpars 
        if (verbose$most) {
          vranfix <- unlist(object$fit$corrPars[c("nu","rho")])
          cat(paste(paste(paste(names(vranfix),"=",signif(vranfix)),collapse=", ")," (estimated by ",itmethod,")\n",sep=""))
        }
      }
      currsurf <- NULL ## garbage management
      object$latestPoints <- nrow(surfaceData)+1-seq_len(nrow(newlogLs)) ## for plots
    }    
    # maximization, new CIs, new MSEs
    object <- MSL(object, CIs=useCI, level=level, verbose=verbose$most)
    if (inherits(object,"SLik_j")) {attr(object$logLs,"freq_good") <- freq_good} 
    MSEs <- object$RMSEs^2
    if (verbose$movie) plot(object,...)
    it <- it+1
    newsimuls <- NULL ## essential for the test
  } ## end while() loop
  if (verbose$most) {
    plot(object,...)
    if  ( any(is.na(RMSEs))) {
      cat("Target precision could not be ascertained.\n")
    } else if (is.null(object$RMSEs)) {
      cat("Precision has not been evaluated.\n")
    } else if  ( any(RMSEs>precision) ) {cat("\nTarget precision does not appear to be reached.\n")}
  }
  return(object)
}