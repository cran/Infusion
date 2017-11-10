## fast version of dmvnorm that uses precomputed version of chol(sigma) and minimal checking
.fast_dmvnorm <- function (x, mean, solve_t_chol_sigma, log = FALSE) {
  if (is.vector(x)) 
    dim(x) <- c(1,length(x)) # x <- matrix(x, ncol = length(x))
  tmp <- solve_t_chol_sigma %*% (t(x) - mean) #backsolve(chol_sigma, t(x) - mean, transpose = TRUE)
  rss <- colSums(tmp^2)
  logretval <- sum(log(diag(solve_t_chol_sigma))) - 0.5 * ncol(x) * log(2 * pi) - 0.5 * rss
  names(logretval) <- rownames(x)
  if (log) 
    logretval
  else exp(logretval)
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
  if (!is.list(verbose)) verbose <- as.list(verbose)
  if (is.null(names(verbose))) names(verbose) <- c("most","final")[seq_len(length(verbose))]
  if (is.null(verbose$most)) verbose$most <- interactive()
  if (is.null(verbose$final)) verbose$final <- FALSE
  statNames <- names(stat.obs)
  if (length(unique(colnames(data)))<ncol(data)) {
    stop("Some columns of the 'data' have identical names. Please correct.")
  }
  if (is.null(allPars <- attr(data,"allPars"))) {  
    allPars <- setdiff(colnames(data),statNames) ## first iteration
  } ## else do not add the cumul_iter variable to the parameters !
  whichvar <- apply(data[,allPars,drop=FALSE],2L,function(v) length(unique(range(v)))>1L)
  fittedPars <- allPars[whichvar]
  fixedPars <- data[1L,!whichvar,drop=FALSE] ## has both names and values
  rownames(fixedPars) <- NULL## avoids warning on rownames when later used in cbind()
  if (! is.list(nbCluster)) nbCluster <- list(jointdens=nbCluster, pardens=nbCluster)
  nbCluster <- lapply(nbCluster,eval.parent)
  models <- mixmodGaussianModel(listModels=Infusion.getOption("mixmodGaussianModel"))
  if (verbose$most) cat(paste("Densities modeling:",nrow(data),"points;\n"))
  locarglist <- list(data=data[,c(fittedPars,statNames)],nbCluster=nbCluster$jointdens, 
                     strategy = mixmodStrategy(seed=123), models=models)
  jointdens <- try(do.call("mixmodCluster",locarglist),silent = TRUE)
  if (inherits(jointdens,"try-error")) {
    locarglist$nbCluster <- eval(Infusion.getOption("nbCluster"))
    jointdens <- do.call("mixmodCluster",locarglist)
  } 
  if (length(nbCluster)==1L) {
    jointdens <- jointdens@results[[1L]]
  } else { jointdens <- .get_best_clu_by_AIC(jointdens) }
  if (verbose$most) cat(paste(" joint density modeling:",jointdens@nbCluster,"clusters;\n"))
  locarglist <- list(data=data[,fittedPars,drop=FALSE], nbCluster=nbCluster$pardens, 
                     strategy = mixmodStrategy(seed=123), models=models)
  pardens <- try(do.call("mixmodCluster",locarglist),silent = TRUE)
  if (inherits(pardens,"try-error")) {
    locarglist$nbCluster <- eval(Infusion.getOption("nbCluster"))
    pardens <- do.call("mixmodCluster",locarglist)
  }                               
  if (length(nbCluster)==1L) {
    pardens <- pardens@results[[1L]]
  } else { pardens <- .get_best_clu_by_AIC(pardens) }
  if (verbose$most) cat(paste(" parameter modeling:",pardens@nbCluster,"clusters.\n"))
  #plotCluster(pardens,data=data[,fittedPars]) to plot @results[[1L]]
  # Rmixmod::plot(<mixmodCluster object>) ith Rmixmod:: to plot from any envir, not only the global one
  resu <- list(jointdens=jointdens,pardens=pardens)
  resu$solve_t_chol_sigma_lists <- list(
    pardens=lapply(resu$pardens@parameters["variance"], function(mat) solve(t(chol(mat)))),
    jointdens=lapply(resu$jointdens@parameters["variance"], function(mat) solve(t(chol(mat))))
  ) 
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
  resu$LOWER <- apply(data[,fittedPars,drop=FALSE],2L,min) ## potentially from an external source such as same-name attr from the input object 
  resu$UPPER <- apply(data[,fittedPars,drop=FALSE],2L,max) # ... used by confintAll
  resu$lower <- apply(data[,fittedPars,drop=FALSE],2L,min) ## really this
  if (any(is.na(c(resu$lower,resu$upper)))) stop("NA(s) in c(lower,upper))")
  resu$upper <- apply(data[,fittedPars,drop=FALSE],2L,max) # ... used by MSL -> optim
  attr(resu,"Infusion.version") <- packageVersion("Infusion")
  class(resu) <- c("SLik_j",class(resu))
  predEDF <- predict(resu,newdata=data[,fittedPars,drop=FALSE])
  wm <- which.max(predEDF)
  resu$optrEDF <- list(par=data[wm,fittedPars,drop=FALSE],value=predEDF[wm])
  return(resu)
}


## infer_SLik_joint -> predict -> predict.SLik_j -> predict.MixmodResults
# FR->FR Il faut que je relise le code de densityMixmod, wrapper pour mixmodCluster
# qui gere les proba masses, et peuxt être pertinent pour joint densities...
predict.MixmodResults <- function (object, newdata,log=TRUE, solve_t_chol_sigma_list,...) {
  Sobs_activeBoundaries <- atb <- freqs <- NULL ## FR->FR tempo
  nbCluster <- object@nbCluster
  if (is.null(nrow(newdata)) ) newdata <- t(as.matrix(newdata)) # as.matrix keeps names
  if (log) { 
    if (nbCluster>0L) {
      density <- matrix(nrow=nrow(newdata),ncol=nbCluster)
      for (k in 1:nbCluster) {
        density[,k] <- log(object@parameters["proportions", k]) + 
          .fast_dmvnorm(newdata, object@parameters["mean", k], solve_t_chol_sigma= solve_t_chol_sigma_list[[k]],log=log)
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
          .fast_dmvnorm(newdata, object@parameters["mean", k], solve_t_chol_sigma_list[[k]],log=log)
      }
      mixture <- rowSums(density) ## sum(Li) 
      if ( !is.null(Sobs_activeBoundaries)) mixture <- mixture*freqs
    } else mixture <- freqs
  }
  return(mixture)
}

predict.SLik_j <- function (object, 
                            newdata, ## requests new fittedPars values! 
                            # newdata=object$logLs[,object$colTypes$fittedPars,drop=FALSE], 
                            log=TRUE, 
                            which="condvaldens",
                            ...) {
  if (is.null(nrow(newdata)) ) newdata <- t(as.matrix(newdata)) # as.matrix keeps names
  if (is.null(colnames(newdata))) colnames(newdata) <- object$colTypes$fittedPars
  tstat.obs <- t(attr(object$logLs,"stat.obs")) ## FR->FR rethink this for all classes of objects
  # not useful bc predict.MixmodResults does not handle it:
  #newdata <- cbind(newdata, t(replicate(nrow(newdata),attr(object$logLs,"stat.obs"))))
  pointpredict <- function(point)  {
    if (which!="parvaldens") jointvaldens <- predict(object$jointdens,
                                                     newdata=cbind(t(point),tstat.obs),
                                                     solve_t_chol_sigma_list=object$solve_t_chol_sigma_lists$jointdens,
                                                     log=log,...)
    if (which=="jointvaldens") return(jointvaldens)
    parvaldens <- predict(object$pardens,newdata=point,
                          solve_t_chol_sigma_list=object$solve_t_chol_sigma_lists$pardens,
                          log=log,...) 
    if (which=="parvaldens") return(parvaldens)
    if (log) {
      condvaldens <- jointvaldens - parvaldens
    } else condvaldens <- jointvaldens/parvaldens
    return(condvaldens)
  }
  apply(newdata,1L,pointpredict)
}

confint.SLik_j <- function(object, parm, ## parm is the parameter which CI is sought 
                         level=0.95, verbose=interactive(),fixed=NULL,which=c(TRUE,TRUE),...) {
  .confintAll(object=object, parm=parm, ## parm is the parameter which CI is sought 
             givenmax = object$MSL$maxlogL,
             level= - qchisq(level,df=1)/2, ## df=1 for 1D profile; /2 to compare to logLR rather than 2 logLR  
             verbose=verbose,fixed=fixed,which=which,...)
}

refine.SLik_j <- function(object,method=NULL,...) {
  if (is.null(method)) method <- "mixmodCluster"
  refine.default(object, surfaceData=object$logLs, method=method, ...)
}

print.SLik_j <-function(x,...) {summary.SLik_j(x,...)} 

.rparam_from_SLik_j <- function(object, size=NULL,fittedPars=NULL,level) {
  if (is.null(fittedPars)) fittedPars <- object$colTypes$fittedPars
  CI_LRstat <- qchisq(level,df=length(fittedPars)) 
  size <- attr(object$logLs,"n_first_iter")
  if (is.null(size)) size <- length(which(object$logLs$cumul_iter==1)) ## FR->FR temporary back compat
  if ( ! is.null(freq_good <- attr(object$logLs,"freq_good"))) size <- ceiling(size/freq_good)
  trypoints <- .simulate.MixmodResults(object$jointdens, nsim=1L, size=size) ## simulates joint
  trypoints <- trypoints[,seq_len(length(fittedPars)),drop=FALSE] ##  marginal=current distr of param
  colnames(trypoints) <- fittedPars
  prior <- predict(object,trypoints,which="parvaldens") ## slow
  logLs <- predict(object,trypoints) ## slow
  criterion <- logLs-prior
  threshold <- max(criterion)-CI_LRstat ## 
  regul <- pmax(threshold-CI_LRstat, ## sets a minimal weight=exp(-CI_LRstat)
                pmin(threshold,criterion)) ## sets a maximal weight = exp(0)
  weights <- exp(regul-threshold) 
  if (FALSE) { ## ring sampling. yields poor clustering.
    locsd <- 1/5
    weights <- weights*sapply(criterion,function(v) max(dnorm(v-c(0,CI_LRstat/2),sd = locsd)))/(0.3989422804014327*locsd)
  }
  good <- (runif(n=length(weights))<weights) ## samples uniformly in the current fitted distrib of params with logL>threshold
  freq_good <- length(which(good))/size
  if (freq_good < 1/5) {
    good <- order(weights,decreasing=TRUE)[seq(size/5)]
    freq_good <- 1/5 ## to avoid explosion of # of sampled points
  }
  trypoints <- trypoints[good,,drop=FALSE]
  trypoints <- (cbind(trypoints,object$colTypes$fixedPars)) ## add fixedPars for simulation
  trypoints <- trypoints[,object$colTypes$allPars,drop=FALSE] ## column reordering
  trypoints <- data.frame(trypoints)
  for (vv in fittedPars) trypoints <- trypoints[trypoints[,vv]>object$LOWER[vv],,drop=FALSE]
  for (vv in fittedPars) trypoints <- trypoints[trypoints[,vv]<object$UPPER[vv],,drop=FALSE]
  return(list(trypoints=trypoints, freq_good=freq_good))
}




## FR->FR define rparam as generic with methods SLik_j, etc ?
.bad_rparam_from_SLik_j <- function(object, size=1L,CIweight=Infusion.getOption("CIweight")) {
  locpredict <- function(x) {predict.SLik_j(object,newdata=x)}
  MSLE <- object$MSL$MSLE
  fittedNames <- names(MSLE)
  RMSEs <- object$RMSEs
  #
  nvec <- .make_n(RMSEs=RMSEs, fittedPars=fittedNames, n=size, CIweight=CIweight)
  #
  trypoints <- matrix(NA,ncol=length(MSLE),nrow=0L)
  if (nvec["MSL"]>0L) {
    hess <- hessian(locpredict,x=MSLE)
    sigma <- - solve(hess) # *(RMSEs["MSL"])
    rmv <- rmvnorm(nvec["MSL"],mean=MSLE,sigma=sigma) # 1st sampling
    trypoints <- rbind(trypoints,rmv) ## centers points on the MSLE 
  }
  CIpoints <- object$CIobject$bounds
  for (kt in seq_len(NROW(CIpoints))) {
    focalpt <- CIpoints[kt,] 
    locst <- rownames(CIpoints)[kt]
    if (nvec[locst]>0L) {
      dotpos <- regexpr(".",locst,fixed=TRUE)[1L]
      parm <- substring(locst,dotpos+1L)
      parpos <- which(fittedNames==parm)
      curv <- hessian(locpredict,x=focalpt)
      curv[parpos,] <- 0
      curv[,parpos] <- 0
      curv[parpos,parpos] <- -grad(locpredict,x=focalpt)[parpos]^2
      sigma <- - solve(curv) # *(RMSEs[locst])
      rmv <- rmvnorm(nvec[locst],mean=focalpt,sigma=sigma)
      trypoints <- rbind(trypoints,rmv)
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
        focalpt[parm] <- bnd*0.999+0.001*focalpt[parm] ## create refpoint near boundary 
        parpos <- which(fittedNames==parm)
        curv <- hessian(locpredict,x=focalpt)
        curv[parpos,] <- 0
        curv[,parpos] <- 0
        curv[parpos,parpos] <- -grad(locpredict,x=focalpt)[parpos]^2
        sigma <- - solve(curv) ## proceeds as if RMSE = 1
        rmv <- rmvnorm(nvec[st],mean=focalpt,sigma=sigma)
        trypoints <- rbind(trypoints,rmv)
      }
    }
  }
  colnames(trypoints) <- fittedNames
  for (vv in fittedNames) trypoints <- trypoints[trypoints[,vv]>object$LOWER[vv],]
  for (vv in fittedNames) trypoints <- trypoints[trypoints[,vv]<object$UPPER[vv],]
  trypoints <- (cbind(trypoints,object$colTypes$fixedPars)) ## add fixedPars for simulation
  trypoints <- trypoints[,object$colTypes$allPars,drop=FALSE] ## column reordering
  return(trypoints)
}
#debug(rparam_from_SLik_j)
#rparam_from_SLik_j(mjd,7)

profile.SLik_j <- function(fitted,...) profile.SLik(fitted=fitted,...) ## no need for two distinct methods here


.boot.SLik_j <- function(object,nsim=2L,verbose = TRUE) {
  if (nsim<2L) stop("'slik' must be > 1")
  bootrepls <- .simulate.MixmodResults(object$jointdens, nsim=nsim, size=nrow(object$logLs))
  boo <- object$logLs[,with(object$colTypes,c(fittedPars,statNames))]
  it <- 0L
  prevmsglength <- 0L
  boo_SLik_joint <- function(simul) {
    it <<- it+1
    locmess <- paste("boostrap replicate",it,"of",nsim,"   ")
    if (verbose) {prevmsglength <<- .overcat(locmess,prevmsglength = prevmsglength)}
    boo[] <- simul
    densv <- infer_SLik_joint(boo,stat.obs=attr(object$logLs,"stat.obs"),
                              nbCluster=list(jointdens=object$jointdens@nbCluster,
                                             pardens=object$pardens@nbCluster),
                              verbose=list(most=FALSE,final=FALSE))
    return(densv)
  }
  resu <- lapply(bootrepls,boo_SLik_joint)
  invisible(resu)
}

.RMSEwrapper.SLik_j <- function(object, CIpoints=object$CIobject$bounds, useCI=TRUE, nsim=10L,verbose=interactive()) {
  bootrepls <- .boot.SLik_j(object,nsim=nsim,verbose=verbose)
  if( useCI && ! is.null(CIpoints) ) {
    locdata <- data.frame(rbind(MSLE=object$MSL$MSLE,CIpoints))
    covmat <- cov( do.call(rbind,lapply(bootrepls,predict,newdata=locdata))) ## covmat of predict=of logL
    MSEs <- c(MSL=covmat[1,1],diag(covmat[-1,-1,drop=FALSE])+covmat[1,1]-2*covmat[1,-1])
  } else {
    locdata <- object$MSL$MSLE
    MSEs <- structure( var( unlist(lapply(bootrepls,predict,newdata=locdata)))  ,names="MSL")
  }
  RMSEs <- sqrt(MSEs)
  if (length(MSEs)>1L) {
    headline <-  "* RMSEs of summary log-L maximum and of its ratio at CI bounds: *\n"
  } else {
    headline <-  "* RMSEs of summary log-L maximum: *\n"
  }
  if(verbose) {
    cat(headline)
    print(RMSEs)
  }
  return(RMSEs)
}

# overcome problems with plot Rmixmod which may need more than the info i na best result
# chercher 
# setMethod(
#   f="plot",
#   signature=c("Mixmod"),
# dans les sources pour savoir ce que Rmixmod fait.
# cf also plotCluster(lik_j$jointdens,data=slik_j$logLs)
.plotClusters <- function(object,which,...) {
  plotCluster(object[[which]],data=object$logLs,...)
}