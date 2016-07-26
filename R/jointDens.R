get_best_clu_by_AIC <- function(cluObject) {
  BICs <- unlist(lapply(cluObject@results,slot,name="criterionValue"))
  logLs <- unlist(lapply(cluObject@results,slot,name="likelihood"))
  dfs <- (2*logLs+BICs)/(log(cluObject@nbSample))
  AICs <- -2*logLs+2*dfs
  return(cluObject@results[[which.min(AICs)]])
}


infer_SLik_joint <- function(data, ## reference table ~ abc
                            stat.obs,
                            logLname=Infusion.getOption("logLname"), ## not immed useful
                            Simulate=attr(data,"Simulate"), ## may be NULL
                            nbCluster=Infusion.getOption("nbCluster"),
                            verbose=list(most=interactive(), ## must be explicitly set to FALSE in knitr examples
                                         final=FALSE)  
) {
  if ( ! is.null( cn <- colnames(stat.obs))) {
    message("Note: 'stat.obs' should be a numeric vector, not a matrix or data.frame.")
    names(stat.obs) <- cn ## minimal patch so that names() can be used, not colnames()
  }
  if (!is.list(verbose)) verbose <-as.list(verbose)
  if (is.null(names(verbose))) names(verbose) <- c("most","final")[seq_len(length(verbose))]
  if (is.null(verbose$most)) verbose$most <- interactive()
  if (is.null(verbose$final)) verbose$final <- FALSE
  statNames <- names(stat.obs)
  if (length(unique(colnames(data)))<ncol(data)) {
    stop("Some columns of the 'data' have identical names. Please correct.")
  }
  allPars <- setdiff(colnames(data),statNames)
  whichvar <- apply(data[,allPars,drop=FALSE],2L,function(v) length(unique(range(v)))>1L)
  fittedPars <- allPars[whichvar]
  fixedPars <- data[1L,!whichvar,drop=FALSE] ## has both names and values
  rownames(fixedPars) <- NULL## avoids warning on rownames when later used in cbind()
  nbCluster <- eval(nbCluster)
  models <- mixmodGaussianModel(listModels=Infusion.getOption("mixmodGaussianModel"))
  if (verbose$most) cat(paste("Densities modeling:",nrow(data),"points;\n"))
  jointdens <- mixmodCluster(data[,c(fittedPars,statNames)],nbCluster=nbCluster, models=models)
  if (length(nbCluster)==1L) {
    jointdens <- jointdens@results[[1L]]
  } else { jointdens <- get_best_clu_by_AIC(jointdens) }
  if (verbose$most) cat(paste(" joint density modeling:",jointdens@nbCluster,"clusters;\n"))
  pardens <- mixmodCluster(data[,fittedPars],nbCluster=nbCluster, models=models)
  if (length(nbCluster)==1L) {
    pardens <- pardens@results[[1L]]
  } else { pardens <- get_best_clu_by_AIC(pardens) }
  if (verbose$most) cat(paste(" parameter modeling:",pardens@nbCluster,"clusters.\n"))
  #plotCluster(pardens,data=data[,fittedPars]) to plot @results[[1L]]
  # Rmixmod::plot(<mixmodCluster object>) ith Rmixmod:: to plot from any envir, not only the global one
  resu <- list(jointdens=jointdens,pardens=pardens)
  attr(resu,"EDFstat") <- "[see this string in infer_SLik_joint]" ## 
  resu$logLs <- structure(data,stat.obs=stat.obs,Simulate=Simulate) ## as in infer_surface.logLs
  # attr(resu,"callfn") <- "infer_SLik_joint"
  resu$projectors <- attr(data,"projectors")
  resu$`Infusion.version` <- packageVersion("Infusion")
  resu$colTypes <- list(allPars=allPars, ## keeps the order of the columns
                                fittedPars=fittedPars,
                                fixedPars=fixedPars,
                                logLname=logLname,
                                statNames=statNames)
  resu$LOWER <- apply(data[,fittedPars],2L,min) ## potentially from an external source such as same-name attr from the input object 
  resu$UPPER <- apply(data[,fittedPars],2L,max) # ... used by confintAll
  resu$lower <- apply(data[,fittedPars],2L,min) ## really this
  if (any(is.na(c(resu$lower,resu$upper)))) stop("NA(s) in c(lower,upper))")
  resu$upper <- apply(data[,fittedPars],2L,max) # ... used by MSL -> optim
  attr(resu,"Infusion.version") <- packageVersion("Infusion")
  class(resu) <- c("SLik_j",class(resu))
  predEDF <- predict(resu,newdata=data[,fittedPars])
  wm <- which.max(predEDF)
  resu$optrEDF <- list(par=data[wm,fittedPars],value=predEDF[wm])
  return(resu)
}


## infer_SLik_joint -> predict -> predict.SLik_j -> predict.MixmodResults
# FR->FR Il faut que je relise le code de densityMixmod, wrapper pour mixmodCluster
# qui gere les proba masses, et peuxt être pertinent pour joint densities...
predict.MixmodResults <- function (object, newdata,log=TRUE, ...) {
  Sobs_activeBoundaries <- atb <- freqs <- NULL ## FR->FR tempo
  nbCluster <- object@nbCluster
  if (is.null(nrow(newdata)) ) newdata <- t(as.matrix(newdata)) # as.matrix keeps names
  if (log) { 
    if (nbCluster>0L) {
      density <- matrix(nrow=nrow(newdata),ncol=nbCluster)
      for (k in 1:nbCluster) {
        density[,k] <- log(object@parameters["proportions", k]) + 
          dmvnorm(newdata, object@parameters["mean", k], sigma= object@parameters["variance",k],log=log)
      }
      maxlogs <- apply(density,1,max)
      normalizedlogs <- apply(density,2L,`-`,maxlogs) ## highest value per row is 0
      ## apply return mess again: 
      if (is.vector(normalizedlogs)) normalizedlogs <- matrix(normalizedlogs,nrow=1)
      mixture <- rowSums(exp(normalizedlogs)) ## exp(normalizedlogs) in (0,1); sum(exp(logLi-maxlog))= exp(-maxlog) sum(exp(logLi))= exp(-maxlog) sum(Li)
      mixture <- log(mixture) + maxlogs ## log(sum(Li))
      if ( !is.null(Sobs_activeBoundaries)) mixture <- mixture + atb*log(freqs)
    } else mixture <- atb*log(freqs)
  } else {
    if (nbCluster>0L) {
      density <- matrix(nrow=nrow(newdata),ncol=nbCluster)
      for (k in 1:nbCluster) {
        density[,k] <- object@parameters["proportions", k] * 
          dmvnorm(newdata, object@parameters["mean", k], sigma= object@parameters["variance",k],log=log)
      }
      mixture <- rowSums(density) ## sum(Li) 
      if ( !is.null(Sobs_activeBoundaries)) mixture <- mixture*freqs
    } else mixture <- freqs
  }
  return(mixture)
}

predict.SLik_j <- function (object, 
                               newdata, ## new fittedPars ! 
                               log=TRUE, ...) {
  if (is.null(nrow(newdata)) ) newdata <- t(as.matrix(newdata)) # as.matrix keeps names
  if (is.null(colnames(newdata))) colnames(newdata) <- object$colTypes$fittedPars
  tstat.obs <- t(attr(object$logLs,"stat.obs")) ## FR->FR rethink this for all classes of objects
  jointvaldens <- predict(object$jointdens,newdata=cbind(newdata,tstat.obs),log=log,...) 
  parvaldens <- predict(object$pardens,newdata=newdata,log=log,...) 
  if (log) {
    condvaldens <- jointvaldens - parvaldens
  } else condvaldens <- jointvaldens/parvaldens
  return(condvaldens)
}

confint.SLik_j <- function(object, parm, ## parm is the parameter which CI is sought 
                         level=0.95, verbose=interactive(),fixed=NULL,which=c(TRUE,TRUE),...) {
  confintAll(object=object, parm=parm, ## parm is the parameter which CI is sought 
             givenmax = object$MSL$maxlogL,
             level= - qchisq(level,df=1)/2, ## df=1 for 1D profile; /2 to compare to logLR rather than 2 logLR  
             verbose=verbose,fixed=fixed,which=which,...)
}

refine.SLik_j <- function(object,...) {
  refine.default(object,surfaceData=object$logLs,...)
}

## FR->FR define rparam as generic with methods SLik_j, etc ?
rparam_from_SLik_j <- function(object,
                                  n=100, ## per focal point
                                  spread=list(MLE=c(1,1/4),CI=c(1,1/4)) ## Sig^2 factor. Would ideally diminish when there is more certainty in focal points 
                                  ) {
  locpredict <- function(x) {predict.SLik_j(object,newdata=x)}
  MSLE <- object$MSL$MSLE
  fittedNames <- names(MSLE)
  hess <- hessian(locpredict,x=MSLE)
  rmv <- rmvnorm(n,sigma=-spread$MLE[1]*solve(hess))
  rmv <- rbind(rmv,rmvnorm(n,sigma=-spread$MLE[2]*solve(hess)))
  trypoints <- sweep(rmv,2L,MSLE,`+`)
  CIpoints <- object$CIobject$bounds
  for (kt in seq_len(NROW(CIpoints))) {
    focalpt <- CIpoints[kt,] 
    locst <- rownames(CIpoints)[kt]
    dotpos <- regexpr(".",locst,fixed=TRUE)[1L]
    parm <- substring(locst,dotpos+1L)
    parpos <- which(fittedNames==parm)
    curv <- hessian(locpredict,x=focalpt)
    curv[parpos,] <- 0
    curv[,parpos] <- 0
    curv[parpos,parpos] <- -grad(locpredict,x=focalpt)[parpos]^2
    rmv <- rmvnorm(n,sigma=-spread$CI[1]*solve(curv))
    rmv <- rbind(rmv,rmvnorm(n,sigma=-spread$CI[2]*solve(curv)))
    trypoints <- rbind(trypoints,sweep(rmv,2L,focalpt,`+`))
  }
  colnames(trypoints) <- fittedNames
  trypoints <- (cbind(trypoints,object$colTypes$fixedPars)) ## add fixedPars for simulation
  trypoints <- trypoints[,object$colTypes$allPars,drop=FALSE] ## column reordering
  return(trypoints)
}
#debug(rparam_from_SLik_j)
#rparam_from_SLik_j(mjd,7)

plot.SLik_j <-function(x,y,filled=FALSE,
                     decorations=NULL, ## ADD to the default decorations
                     ...) {
  object <- x
  fittedPars <- object$colTypes$fittedPars
  maxlogL <- object$MSL$maxlogL
  if (is.null(maxlogL)) stop("plot.SLik plots likelihood ratios, hence it requires the ML to have been computed by MSL(.)\n")
  Ztransf <- function(Z) {exp(Z-maxlogL)}
  np <- length(fittedPars)
  if (np==1L) {
    x <- object$logLs[,fittedPars]
    y <- Ztransf(predict(object,newdata=x)) ## <>1
    plot(x,y,main="Summary Likelihood Ratio",xlab=fittedPars,ylab="SL ratio",
         xlim=range(x),
         ylim=c(0,max(1,y)))
    points(object$MSL$MSLE,1,pch="+")
    if (!is.null(object$CIobject)) {
      yci <- rep(exp(-qchisq(0.95,1)/2),NROW(object$CIobject$bounds))
      points(object$CIobject$bounds,y=yci,pch="+")
    }
  } else if (np==2L) {
    decos <- quote({if (!is.null(object$latestPoints)) {points(object$logLs[object$latestPoints,fittedPars],pch=".",cex=2)};
      if (!is.null(object$CIobject)) {points(object$CIobject$bounds,pch="+",cex=1.5,col="red")};
      points(t(object$MSL$MSLE),pch="+",cex=2);
      decorations}) ## language object
    xy <- object$logLs[,fittedPars]
    spaMMplot2D(x=xy[,1],y=xy[,2],z=Ztransf(predict(object,newdata=xy)),
                color.palette=function(n){spaMM.colors(50,redshift=3)},nlevels=50,
                plot.title=quote(title(main="Summary likelihood ratio surface",
                                       xlab=fittedPars[1],ylab=fittedPars[2])),
                decorations=decos,
                ...)
  } else if (np>2L) {
    intsqrt <- floor(sqrt(np))
    if (intsqrt>1) {loccex.axis <- par("cex.axis")*0.6} else {loccex.axis <- par("cex.axis")}
    ## mfrow marche pas avec rstudio (?) cf OKsmooth::provideDevice
    dev <- getOption("device")
    rstudioMess <-  ( (class(dev)=="character" && dev == "RStudioGD") )  
    if (! rstudioMess) opar <- par(cex.axis=loccex.axis, no.readonly = TRUE)
    #  # if (! rstudioMess) opar <- par(mfrow=c(ceiling(np/intsqrt), intsqrt), cex.axis=loccex.axis, no.readonly = TRUE)
    #  ## cf blackbox::gridfn, makeplot, etc
    # mais en fait migraine necase pas plusieurs filled plots sur un device; ici ca ferait planter car plot.new() -> figure margins too large
    margefrac <- 0.001
    gridsteps <- 40
    np <- length(fittedPars)
    grillelist <- list()
    ranges <- apply(object$logLs[,fittedPars],2,range)
    for (it in seq_len(np-1)) {
      xvar <- fittedPars[it]
      range <- ranges[,xvar]
      lob <- range[1]
      upb <- range[2]
      marge <-  (upb-lob)/(4*gridsteps) ## just enough to see the maximum on the edge
      grillelist[xvar] <- list(seq(lob-marge, upb+marge, length.out=gridsteps)) ## uniform on log scale if relevant
      for (jt in (it+1):np) {
        yvar <- fittedPars[jt]
        fixedPars <- setdiff(fittedPars,c(xvar,yvar))
        grillelist[fixedPars] <- NULL
        fixedVals <- object$MSL$MSLE[fixedPars]
        range <- ranges[,yvar]
        lob <- range[1]
        upb <- range[2]
        marge <-  (upb-lob)/(4*gridsteps) 
        grillelist[yvar] <- list(seq(lob-marge, upb+marge, length.out=gridsteps)) ## uniform on log scale if relevant
        ## oddly, order in grillelist is not always well controlled at this point hence
        grillelist <- grillelist[c(xvar,yvar)] ## simply reorder elements according to this order
        grille <- expand.grid(grillelist) 
        grille <- cbind(grille,t(fixedVals))
        grille <- grille[,fittedPars] ## simply reorder grille elements according to fittedNames order
        z <- predict(object, grille)
        xyz <- as.surface(grillelist, z, order.variables = "xy")
        main <- paste("logL slice for",paste(fixedPars,signif(fixedVals,4),sep="=",collapse=", "))
        varVals <- object$MSL$MSLE[c(xvar,yvar)]
        if (interactive() && ! rstudioMess) plot.new() # provideDevice(bbdefaultPars=TRUE) does not work by itself, cf providePlotFile
        spaMM.filled.contour(xyz$x, xyz$y, xyz$z,xlab=xyz$xlab,ylab=xyz$ylab,main=main,
                             color.palette= function(v) {spaMM.colors(v, redshift = 1/2)},
                             plot.axes={
                               axis(1); axis(2); 
                               points(varVals[xyz$xlab],varVals[xyz$ylab],pch="+",cex=1.5) 
                             }
        )
      } 
    }
    if ( ! rstudioMess) par(opar)
  }
  invisible(object)
}

profile.SLik_j <- function(fitted,...) profile.SLik(fitted=fitted,...) ## no need for two distinct methods here