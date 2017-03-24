## auxiliary fn for sample_volume, tiself for rparam
`sample_vertices` <- function(n,tryn,vertices,object,fixed=NULL,fittedPars,useEI)  {
  if (length(fittedPars)>length(fixed)) { ## profile CIs, or 
    if (length(fixed)>0L) {
      fixedPars <- names(fixed)
      if (is.null(fixedPars)) {stop("'fixed' argument of 'sample_vertices' has no name. Check call.")}
      subHull <- subHullWrapper(vertices=vertices,equality=fixed)
      if (is.null(subHull$vertices)) {
        return(NULL)
      } else {
        vT <- volTriangulationWrapper(as.matrix(subHull$vertices))   
        if (inherits(vT,"try-error")) {
          if (substr(vT[1L],9,20)=="nrow <= ncol") {
            locvert <- unique(subHull$vertices)
            nr <- nrow(locvert)
            if (nr==1L) {
              trypoints <- locvert ## but then there are fewer than tryn of them... 
            } else {
              locfn <- function() {
                locweights <- rsimplex(simplex = rbind(1/nr,diag(nr)))
                colSums(sweep(locvert,1L,locweights,`*`)) ## sweep(...) multiply rows by weights
              }
              trypoints <- t(replicate(tryn, locfn()))
            } 
          } else return(vT) ## unhandled error in volTriangulationWrapper <=> should not occur  
        } else {
          trypoints <- data.frame(rvolTriangulation(tryn,vT))
        }
      }
      trypoints <- cbind(fixed,trypoints,row.names=NULL) 
      colnames(trypoints)[seq_len(length(fixedPars))] <- fixedPars ## util ?
      trypoints <- trypoints[,fittedPars,drop=FALSE] ## reordering
      Qmax <- NULL ## global Qmax not useful for selecting CI candidates
    } else {
      vT <- volTriangulationWrapper(as.matrix(vertices))
      trypoints <- data.frame(rvolTriangulation(tryn,vT))
      Qmax <- object$Qmax
    }
    if (useEI && nrow(trypoints)>n) { ## nrow(trypoints)>n may be FALSE in the degenerate case where trypoints = unique(subHull$vertices)
      expectedImprovement <- calc_EI(object,trypoints,Qmax=Qmax)  
      trypoints <- cbind(trypoints,EI=expectedImprovement)
      trypoints <- trypoints[order(trypoints[,"EI"],decreasing=TRUE)[seq_len(n)],fittedPars,drop=FALSE]
    } else if (tryn>n) {warning("From sample_vertices(): 'tryn'>'n' is uselessly time-consuming.")}
  } else {
    ## in one parameter model there is no tengant subspace of proiled out params to investigate
    ## using EI on CI variable would tend to select points towards the MSLE
    ## hence EI is not used in this case.
    locrange <- range(vertices[,1])
    if (length(fixed)==1L) locrange <- locrange*0.5+fixed[1L]*0.5 ## more vaguely focused sampling
    trypoints <- runif(n,min=locrange[1],max=locrange[2]) ## FR->FR potential for improvement ? analogous to quadrature points?
    if (any(is.na(trypoints))) stop("NA(s) in trypoints")
    trypoints <- data.frame(trypoints)
    vT <- NULL
  }
  return(trypoints)
}


## auxiliary fn for rparam
## j'ai rajouté les bary dans la sortie. Y penser pour migraineLike
`sample_volume` <- function(object,n=6,useEI,vertices=NULL,
                            dlr=NULL, ## NULL => calc.lrthreshold will use a dfault value
                            verbose=interactive(),fixed=NULL,
                            tryn = 30*n ) {
  fittedPars <- object$colTypes$fittedPars
  if (n==0) {    
    resu <- as.data.frame(matrix(nrow = 0, ncol = length(fittedPars)))
    colnames(resu) <- fittedPars
    return(resu)
  } ## ELSE 
  if (useEI) {
    ## so that oldXnew matrices will contain less than spaMM.getOption("ff_threshold")~1e7 elements: 
    maxn <- floor(1e7/nrow(object$fit$data))
    if (maxn <= n) {
      locmess <- paste("From sample_volume(): 'maxn': ",maxn,"<=",n," ('n'). 'n' reduced to")
      n <- ceiling(maxn/10)
      message(paste(locmess,n))
    }
    if (tryn > maxn) {
      locmess <- paste("From sample_volume(): 'tryn' reduced from",tryn,"to",maxn)
      message(locmess)
      tryn <- maxn
    }
  } else tryn <- n  
  if ( ! identical(-Inf,dlr)) {
    lrthreshold <- calc.lrthreshold(object,dlr=dlr,verbose=verbose)  ## with different methods for SLik and SLikp    
    if ( is.null(vertices)) {
      vertices <- object$fit$data[,fittedPars,drop=FALSE]
      pred <- fitted(object$fit)
    } else pred <- predict(object,newdata=vertices)
    uppervertices <- vertices[pred>lrthreshold,,drop=FALSE]
    ## if there is only one uppervertice, sample_vertices will typically sample only one point constructed from fixed and that vertice
    ## => test to avoid that; but more general pb must be caught by the test:
    if (nrow(uppervertices)>ncol(vertices)) { 
      trypoints <- sample_vertices(n,tryn,vertices=uppervertices,object=object,
                                   fixed=fixed,fittedPars=fittedPars,useEI=useEI) ## can be NULL
    } else  trypoints <- NULL
  } else {
    if ( is.null(vertices)) vertices <- object$fit$data[,fittedPars,drop=FALSE]
    trypoints <- NULL
  }
  if (is.null(trypoints)) trypoints <- sample_vertices(n,tryn,vertices=vertices,object=object,
                                                       fixed=fixed,fittedPars=fittedPars,useEI=useEI)
  oldn <- NROW(trypoints)
  trypoints <- unique(trypoints)
  if (NROW(trypoints) < oldn) {warning("From sample_volume(): Suspect event: replicate points generated by sample_vertices(). ")}
  #  attr(trypoints,"vT") <- vT
  colnames(trypoints) <- fittedPars
  return(trypoints) ## fittedPars parameters
}

.make_n <- function(RMSEs, fittedPars, n, CIweight) {
  if (is.null(RMSEs)) RMSEs <- 1e10
  np <- length(fittedPars)
  if(is.null(n)) n <- ( 2+1*(np>1) ) * (2*np+1) # 1st factor => 1 ou 3 profile points for each bound; 2nd => cf max length of MSEs
  CIweight <- CIweight*(n+1)/(n-2)
  wMSEs <- rep(CIweight,(2*np+1))
  names(wMSEs) <- c("MSL",as.vector( rbind(paste("low",fittedPars,sep="."),paste("up",fittedPars,sep="."))))
  wMSEs["MSL"] <- 1 ## relative to CIweight: change CIweight rather than this.
  if ( ! any(is.na(RMSEs))) wMSEs[names(RMSEs)] <- wMSEs[names(RMSEs)] * RMSEs^2 ## weighted MSEs giving double weight to CI bounds
  En <- n*wMSEs/sum(wMSEs)
  nvec <- floor(En)
  En <- En - nvec
  sumEn <- sum(En)
  if (sumEn>0) {
    nsample <- sample(2*np+1,size=n-sum(nvec),prob=En/sumEn,replace=TRUE)
    messy <- table(nsample)
    numnames <- as.numeric(names(messy))
    nvec[numnames] <- nvec[numnames] + messy
  }
  names(nvec) <- names(wMSEs)
  return(nvec)
}

# return value includes fixedPars
rparam <- function(object, n= 1, useEI = list(max=TRUE,profileCI=TRUE,rawCI=FALSE), useCI = TRUE, 
                   verbose = interactive(), tryn=30*n,  
                   level = 0.95, 
                   #auxlevel=1-((1-level)*2/3),
                   CIweight=Infusion.getOption("CIweight")) {
  fittedPars <- object$colTypes$fittedPars  
  nvec <-  .make_n(RMSEs=object$RMSEs, fittedPars=fittedPars, n=n, CIweight=CIweight)
  np <- length(fittedPars)
  #
  pts <- object$fit$data[,object$colTypes$fittedPars,drop=FALSE]
  vT <- volTriangulation(pts)
  fallback <- cbind(object$LOWER,object$UPPER)
  gridinfo <- lapply(object$colTypes$fittedPars,function(par) fallback[par,])
  names(gridinfo) <- object$colTypes$fittedPars
  largerpts <- do.call(expand.grid,gridinfo)
  largervT <- volTriangulation(rbind(pts,largerpts))
  #
  if (useCI) {
    trypoints <- allCIs(object,level=0.1,verbose=FALSE)$bounds ## narrow bracketing of $optr$par by level=0.1
    if (nvec["MSL"]<NROW(trypoints)) { 
      trypoints <- trypoints[sample(nrow(trypoints),nvec["MSL"]),,drop=FALSE]
    } 
  } else trypoints <- NULL
  if (NROW(trypoints)<nvec["MSL"]) {
    whichSimplex <- locatePointinvT(object$MSL$MSLE,vT)
    subvT <- subsimplices.volTriangulation(vT,whichSimplex)
    vertices <- subvT$vertices[unique(as.vector(subvT$simplicesTable)),,drop=FALSE]
    morepoints <- sample_volume(object,n=nvec[1L],useEI=useEI$max,vertices=vertices,
                               verbose=verbose) ## EI near the maximum point
    trypoints <- rbind(trypoints,morepoints)
    #print(trypoints)
  }
    # points targeting CI bounds
  if (useCI) {
    if (is.null(object$CIobject))  object$CIobject <- allCIs(object,level=level,verbose=FALSE)
    auxCIobject <- allCIs(object,level=1-((1-object$CIobject$level)*2/3),verbose=FALSE)
    CIpoints <- object$CIobject$bounds
    ## those for which the bound is present:
    for (kt in seq_len(NROW(CIpoints))) { 
      focalpt <- CIpoints[kt,] ## loses names from 1-col matrix ! Hence:
      if (ncol(CIpoints)==1L) names(focalpt) <- colnames(CIpoints)
      locst <- rownames(CIpoints)[kt]
      dotpos <- regexpr(".",locst,fixed=TRUE)[1L]
      parm <- substring(locst,dotpos+1L)
      # the problem with using sample_volume is the selection of vertices...
      whichSimplex <- locatePointinvT(focalpt,vT,fallback=FALSE)
      if (length(whichSimplex)==0L) {
        whichSimplex <- locatePointinvT(focalpt,largervT,fallback=TRUE)
        subvT <- subsimplices.volTriangulation(largervT,whichSimplex)
        ## FR->FR and maybe also call rsimplex to sample more widely an interesting unsampled region ?
      } else subvT <- subsimplices.volTriangulation(vT,whichSimplex)
      vertices <-subvT$vertices[unique(as.vector(subvT$simplicesTable)),,drop=FALSE]
      fixed <- focalpt[parm]
      if (np>length(fixed)) {
        loc_useEI <- useEI$profileCI
      } else loc_useEI <- useEI$rawCI
      if (locst %in% rownames(auxCIobject$bounds)) {
        n_onspot <- ceiling(nvec[kt+1L]/2)
      } else n_onspot <- nvec[kt+1L]
      bnds <- sample_volume(object,n=n_onspot,useEI=loc_useEI,vertices=vertices,
                           fixed=fixed,verbose=verbose) 
      trypoints <- rbind(trypoints,bnds)
      #
      n_outspot <- nvec[kt+1L]-n_onspot
      if (n_outspot>0L) {
        auxCIpoints <- auxCIobject$bounds
        focalpt <- auxCIpoints[locst,] ## loses names from 1-col matrix ! Hence:
        if (ncol(auxCIpoints)==1L) names(focalpt) <- colnames(auxCIpoints)
        whichSimplex <- locatePointinvT(focalpt,vT,fallback=FALSE)
        if (length(whichSimplex)==0L) {
          whichSimplex <- locatePointinvT(focalpt,largervT,fallback=TRUE)
          subvT <- subsimplices.volTriangulation(largervT,whichSimplex)
        } else subvT <- subsimplices.volTriangulation(vT,whichSimplex)
        vertices <-subvT$vertices[unique(as.vector(subvT$simplicesTable)),,drop=FALSE]
        fixed <- focalpt[parm]
        bnds <- sample_volume(object,n=n_outspot,useEI=loc_useEI,vertices=vertices,
                             fixed=fixed,verbose=verbose) ## a noter que 
        trypoints <- rbind(trypoints,bnds)
      }
    } 
    # those for which no CI bound is available:
    for (st in object$CIobject$missingBounds) {
      if (nvec[st]>0L) {
        dotpos <- regexpr(".",st,fixed=TRUE)[1L]
        parm <- substring(st,dotpos+1L)
        side <- substring(st,0,dotpos-1L)
        focalpt <- object$MSL$MSLE
        if (side=="low") {
          bnd <- object$LOWER[parm]
        } else bnd <- object$UPPER[parm]
        if (abs(bnd-focalpt[parm])>(object$UPPER[parm]-object$LOWER[parm])/1e4) { ## if MSLE not at boundary
          focalpt[parm] <- bnd*0.999+0.001*focalpt[parm] ## create refpoint near boundary and sample its simplex
          whichSimplex <- locatePointinvT(focalpt,largervT,fallback=TRUE)
          subvT <- subsimplices.volTriangulation(largervT,whichSimplex)
          rpts <- t(replicate(nvec[st], rsimplex(simplex=subvT$vertices[subvT$simplicesTable,])))
          trypoints <- rbind(trypoints,rpts)
        }
      }
    }
  }
  # we have 'n' points
  if (nrow(trypoints)>0L) { ## bc rbind has an aberrant naming behaviour otherwise
    withRepl <- c(seq_len(nrow(trypoints)),sample(nrow(trypoints),1L)) ## original 'n' plus one repliate => 'n+1' 
    trypoints <- trypoints[withRepl,,drop=FALSE]
    trypoints <- rbind(trypoints,object$MSL$MSLE) ## fittedPars only => 'n+2'
  } else trypoints <- t(object$MSL$MSLE) ## (allows n=0)
  if (verbose) {
    cat("New design points:\n")
    print(trypoints)
  }
  ## in extreme cases MSL$MSLE may be constant over iterations. We check all points to retain at most two replicates:
  newold <- proxy::dist(trypoints,object$logLs[,fittedPars,drop=FALSE])
  trypoints <- trypoints[apply(newold,1,function(line) {length(which(line==0))<2L}),,drop=FALSE] 
  trypoints <- cbind(trypoints,object$colTypes$fixedPars) ## add fixedPars for simulation
  trypoints <- trypoints[,object$colTypes$allPars,drop=FALSE] ## column reordering
  return(trypoints) ## typically 'n+2' with 'n+1' unique and one replicate, 
                    ## but possibly '1' (if n=0)or 'n+1' (if excess replicates were removed)
}

