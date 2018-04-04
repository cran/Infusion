.densityMclust <- function(data,stat.obs,nbCluster=Infusion.getOption("nbCluster")) {
  # stat.obs useful for boundaries and to handle degenerate distributions:
  nbCluster <- eval(nbCluster)
  statNames <- colnames(data)
  boundaries <- attr(stat.obs,"boundaries") 
  Sobs_activeBoundaries <- simuls_activeBoundaries <- freq <- NULL
  if ( ! is.null(boundaries) ) {
    simulsMatches <- t(apply(data[,names(boundaries),drop=FALSE],1L,`==`,y=boundaries))
    if (length(boundaries)==1L) { ## apply() fix
      simulsMatches <- t(simulsMatches)
      colnames(simulsMatches) <- names(boundaries)
    }
    simulsAtSomeBoundaries <- apply(simulsMatches,1L,any)
    SobsMatches <- (stat.obs[names(boundaries)]==boundaries)
    SobsAtSomeBoundaries <- any(SobsMatches)
    if (any(SobsAtSomeBoundaries)) {
      Sobs_activeBoundaries <- boundaries[which(SobsAtSomeBoundaries)]
      activeNames <- names(Sobs_activeBoundaries)
      ## keep only simuls that match all the boundaries that stat.obs matches:
      goodrows <- which(apply(simulsMatches[,activeNames,drop=FALSE],1L,all))
      if (length(goodrows)==0L) {
        resu <- structure(resu=new("mixmodCluster",nbCluster=0L),
                          statNames=statNames,
                          Sobs_activeBoundaries=Sobs_activeBoundaries,
                          freq=1/(2*(nrow(data+1))))
        class(resu) <- "dMixmod"
        return(resu)
      } ## ELSE
      ## estimate proba that all (SobsMatches) are observed in the simulation 
      freq <- (length(goodrows)+1/2)/(nrow(data)+1) ## log(freq) is estimator of log(p) with bias ~1/n^2 
      data <- data[goodrows,setdiff(statNames,activeNames),drop=FALSE] ## conditional on match, lower-dimensional
    } else if (any(simulsAtSomeBoundaries)) {
      simuls_activeBoundaries <- boundaries[which(apply(simulsMatches,2L,any))]
      ## keep only simuls that match none of the the boundaries:
      goodrows <- which( ! apply(simulsMatches[,names(simuls_activeBoundaries),drop=FALSE],1L,any)) 
      ## estimate proba that no (simulsMatches) is observed in the simulation 
      freq <- (length(goodrows)+1/2)/(nrow(data)+1) ## log(freq) is estimator of log(p) with bias ~1/n^2 
      data <- data[goodrows,,drop=FALSE] ## conditional on no match, full-dimensional
    }
  } 
  ## 
  checkranges <- apply(data,2,range)
  rangewidths <- checkranges[2L,]-checkranges[1L,]
  degenerates <- rangewidths<1e-6
  if (any(degenerates)) {
    which_dg <- which(degenerates)
    which_out <- (stat.obs[which_dg]>checkranges[2L,which_dg]) | (stat.obs[which_dg]<checkranges[1L,which_dg])
    if (any(which_out)) {
      for (dg in which_dg) {
        dst <- max(abs(stat.obs[dg]-data[,dg]))
        data[,dg] <- data[,dg] + rnorm(n=nrow(data),sd=dst/5) ## makes it reasonably unlikely that obs will be within noised data
      }
    } else  {
      ## stat.obs in degenerate simuls but mixmodCluster needs non-degenerate distribution 
      for (dg in which_dg) data[,dg] <- data[,dg] + rnorm(n=nrow(data),sd=1e-06)
    }
  }
  # 
  #  x <- mixmodCluster(as.data.frame(data),seq(2*log(nrow(data))),strategy = mixmodStrategy(seed=123))
  # statdoc of mixmod recommends ceiling(nrow(data)^0.3) and refers to Bozdogan93
  if (length(nbCluster)==1L) {
    resu <- suppressWarnings(.mixclustWrap("Mclust",
                          list(data=as.data.frame(data),modelNames=Infusion.getOption("mclustModel"), G=nbCluster,verbose=FALSE),
                          pack="mclust"))
    if (is.null(resu)) { ## If mclust fails bc of a singular gaussian component, it returns NULL...
      # try decreasing the number of clusters until mclust works
      for (nb in rev(seq_len(min(nbCluster)-1L))) {
        resu <- suppressWarnings(.mixclustWrap("Mclust",
                                               list(data=as.data.frame(data),modelNames=Infusion.getOption("mclustModel"), G=nb,verbose=FALSE),
                                               pack="mclust"))
        if ( ! is.null(resu) ) break
      }
    } 
  } else { 
    models <- vector("list",length(nbCluster))
    for (it in seq_along(nbCluster)) models[[it]] <- .mixclustWrap("Mclust",
                                                                   list(data=as.data.frame(data),modelNames=Infusion.getOption("mclustModel"), G=nbCluster[it],verbose=FALSE),
                                                                   pack="mclust")
    resu <- .get_best_mclust_by_AIC(models) 
  }
  if (resu$modelName=="VEV") { 
    ## then cov matrices are proportional, wecheck they are not too heterogeneous
    varmats <- resu$parameters$variance$sigma
    vars11 <- apply(varmats,3L,`[`,x=1L,y=1L)
    rangelv <- range(log(vars11))
    if ( rangelv[2L]-rangelv[1L]> 30) {
      message("Inferred Gaussian components have very heterogeneous variances.")
      message("This can occur when some summary statistic takes one particular value with high probability,")
      message("which will be poorly fitted by Mclust().")
      warning("Result is suspect.")
    } 
  } 
  resu <- structure(resu,statNames=statNames,simuls_activeBoundaries=simuls_activeBoundaries,
                    Sobs_activeBoundaries=Sobs_activeBoundaries,freq=freq)
  class(resu) <- "dMclust"
  return(resu)
}

.get_best_mclust_by_AIC <- function(cluObject) {
  BICs <- logLs <- numeric(length(cluObject))
  for (it in seq_along(cluObject)) {
    BICs[it] <- cluObject[[it]]$BIC
    logLs[it] <- cluObject[[it]]$loglik
  }
  dfs <- (2*logLs+BICs)/(log(cluObject[[1]][["n"]]))
  AICs <- -2*logLs+2*dfs
  return(cluObject[[which.min(AICs)]])
}

predict.dMclust <- function(object,
                            newdata, ## should be (reformatted) stat.obs (in canned procedures)
                            tcstat.obs=NULL, ## to avoid checks of arguments
                            log=FALSE,...) {
  statNames <- attr(object,"statNames")
  if (is.null(tcstat.obs)) {
    ns <- length(statNames)
    if (is.null(dim(newdata))) { ## less well controlled case, but useful for maximization (which is not performed in canned procedures)
      newdata <- matrix(newdata,nrow=1L)
      if (ns==ncol(newdata)) {
        colnames(newdata) <- statNames
      } else {
        stop(paste("(!) newdata has incorrect length. It should match the following variables:\n",
                   paste(statNames,collapse=" ")))
      }
    } 
    if ( !is.null(Sobs_activeBoundaries <- attr(object,"Sobs_activeBoundaries"))) {
      boundsdata <-  newdata[,names(Sobs_activeBoundaries),drop=FALSE]
      atb <- apply(boundsdata,1L,`==`,y=Sobs_activeBoundaries)
      if (!all(atb)) {warning("'predict.dMclust' cannot compute joint out-of-boundary density from conditional at-boundary density. ")}
      freq <- attr(object,"freq")
      freqs <- atb*freq+(1-atb)*(1-freq) ## uses 1-freq instead of the unknown density of the variable(s) in boundaries 
      densitydata <- newdata[,setdiff(statNames,names(Sobs_activeBoundaries)),drop=FALSE]
    } else if ( !is.null(simuls_activeBoundaries <- attr(object,"simuls_activeBoundaries"))) {
      ## do not reduce data in this case
      densitydata <- newdata 
      ## only for the warning:
      boundsdata <-  newdata[,names(simuls_activeBoundaries),drop=FALSE]
      atb <- apply(boundsdata,1L,`==`,y=simuls_activeBoundaries)
      if (any(atb)) {
        warning("'predict.dMclust' cannot compute conditional at-boundary density from joint out-of-boundary density.")
        # return value is the same as for predict(,tcstat.obs=<newdata>) 
      }
    } else densitydata <- newdata[,statNames,drop=FALSE]
  } else { 
    if ( !is.null(Sobs_activeBoundaries <- attr(object,"Sobs_activeBoundaries"))) {
      atb <- TRUE
      freqs <- attr(object,"freq")
      densitydata <- tcstat.obs[,setdiff(statNames,names(Sobs_activeBoundaries)),drop=FALSE]
    } else densitydata <- tcstat.obs  ## [,statNames,drop=FALSE]
  } 
  nbCluster <- object$G
  if (log) { 
    if (nbCluster>0L) {
      density <- matrix(nrow=nrow(densitydata),ncol=nbCluster)
      for (k in 1:nbCluster) {
        density[,k] <- object$parameters$pro[k] * 
          dmvnorm(densitydata, object$parameters$mean[,k], sigma= object$parameters$variance$sigma[,,k],log=log)
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
      density <- matrix(nrow=nrow(densitydata),ncol=nbCluster)
      for (k in 1:nbCluster) {
        density[,k] <- object$parameters$pro[k] * 
          dmvnorm(densitydata, object$parameters$mean[,k], sigma= object$parameters$variance$sigma[,,k],log=log)
      }
      mixture <- rowSums(density) ## sum(Li) 
      if ( !is.null(Sobs_activeBoundaries)) mixture <- mixture*freqs
    } else mixture <- freqs
  }
  return(mixture)
}

