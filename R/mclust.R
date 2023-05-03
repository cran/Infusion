.densityMclust <- function(data,stat.obs,nbCluster=seq_nbCluster(nr=nrow(data))) {
  # stat.obs useful for boundaries and to handle degenerate distributions:
  nbCluster <- get_nbCluster_range(projdata=data, nbCluster=nbCluster)
  statNames <- colnames(data)
  ## 
  blob <- .process_boundaries(data, stat.obs, boundaries=attr(stat.obs,"boundaries"))
  if ( ! nrow(blob$data)) { # ie if length(blob$goodrows)==0L in .process_boundaries()
    ## code has long been invalid, but not tested. See further comments on parallel code for .densityMixmod()
    # This still ooks as an ugly patch since Rmixmod is potentially used in this mClust based function.
    resu <- new(structure("dMixmod", package="Infusion")) # line individually OK; no 'data' nor 'nbCluster' slots
    resu@varNames <- blob$varNames
    resu@statNames <- statNames
    resu@simuls_activeBoundaries <- blob$simuls_activeBoundaries
    resu@Sobs_activeBoundaries <- blob$Sobs_activeBoundaries
    resu@freq <- blob$freq # 1/(2*(nrow(data+1)))
  } else {
    # data lines matching the 'boundaries' values have now been removed. There may still be NA's in the data (notably, lines of NA's). The doc explains
    # "the inference method ... must thus ignore all empirical distributions including NA/NaN/inf.": 
    data <- .check_data_post_boundaries(data)
    ## 
    checkranges <- sapply(data,range)
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
      resu <- suppressWarnings(.do_call_wrap("Mclust",
                                             list(data=as.data.frame(data),modelNames=Infusion.getOption("mclustModel"), G=nbCluster,verbose=FALSE),
                                             pack="mclust"))
      if (is.null(resu)) { ## If mclust fails bc of a singular gaussian component, it returns NULL...
        # try decreasing the number of clusters until mclust works
        for (nb in rev(seq_len(min(nbCluster)-1L))) {
          resu <- suppressWarnings(.do_call_wrap("Mclust",
                                                 list(data=as.data.frame(data),modelNames=Infusion.getOption("mclustModel"), G=nb,verbose=FALSE),
                                                 pack="mclust"))
          if ( ! is.null(resu) ) break
        }
      } 
    } else { 
      models <- vector("list",length(nbCluster))
      for (it in seq_along(nbCluster)) models[[it]] <- .do_call_wrap("Mclust",
                                                                     list(data=as.data.frame(data),modelNames=Infusion.getOption("mclustModel"), G=nbCluster[it],verbose=FALSE),
                                                                     pack="mclust")
      resu <- .get_best_mclust_by_IC(models) 
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
    resu <- structure(resu,
                      statNames=statNames,
                      simuls_activeBoundaries=blob$simuls_activeBoundaries,
                      Sobs_activeBoundaries=blob$Sobs_activeBoundaries,
                      freq=blob$freq)
    class(resu) <- c("dMclust", class(resu))
  }
  resu
}

.get_best_mclust_by_IC <- function(cluObject, which=Infusion.getOption("criterion")) {
  if (inherits(cluObject,"try-error")) return(cluObject) ## passes original error info rather than hiding it under another error
  if (length(cluObject)==1L) return(cluObject[[1L]])
  BICs <- logLs <- numeric(length(cluObject))
  for (it in seq_along(cluObject)) {
    BICs[it] <- cluObject[[it]]$BIC
    logLs[it] <- cluObject[[it]]$loglik
  }
  if (which=="BIC") return(cluObject[[which.min(BICs)]])
  # ELSE
  dfs <- (2*logLs+BICs)/(log(cluObject[[1]][["n"]]))
  AICs <- -2*logLs+2*dfs
  return(cluObject[[which.min(AICs)]])
}

predict.dMclust <- function(object,
                            newdata, ## should be (reformatted) stat.obs (in canned procedures)
                            tcstat.obs=NULL, ## to avoid checks of argument
                            solve_t_chol_sigma_list,
                            log=FALSE,...) {
  statNames <- attr(object,"statNames")
  if (is.null(tcstat.obs)) {
    ns <- length(statNames)
    if (is.null(dim(newdata))) { ## less well controlled case, but useful for maximization (which is not performed in canned procedures)
      if ((ns <- length(statNames)) != length(newdata)) {
        stop(paste("(!) newdata has incorrect length. It should match the following variables:\n",
                   paste(statNames,collapse=" ")))
      } else {
        datanames <- names(newdata) # before they are erased by the next dim() assignment (!)
        dim(newdata) <- c(1L, ns)
        if (is.null(datanames)) {
          colnames(newdata) <- statNames
        } else colnames(newdata) <- datanames
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
          .fast_dmvnorm(densitydata, object$parameters$mean[,k], solve_t_chol_sigma= solve_t_chol_sigma_list[[k]],log=log)
      }
      maxlogs <- matrixStats::rowMaxs(density)
      normalizedlogs <- apply(density,2L,`-`,maxlogs) ## highest value per row is 0
      ## apply return mess again: 
      if (is.vector(normalizedlogs))  dim(normalizedlogs) <- c(1L, length(normalizedlogs))
      mixture <- rowSums(exp(normalizedlogs)) ## exp(normalizedlogs) in (0,1); sum(exp(logLi-maxlog))= exp(-maxlog) sum(exp(logLi))= exp(-maxlog) sum(Li)
      mixture <- log(mixture) + maxlogs ## log(sum(Li))
      if ( !is.null(Sobs_activeBoundaries)) mixture <- mixture + atb*log(freqs)
    } else mixture <- atb*log(freqs)
  } else {
    if (nbCluster>0L) {
      density <- matrix(nrow=nrow(densitydata),ncol=nbCluster)
      for (k in 1:nbCluster) {
        density[,k] <- object$parameters$pro[k] * 
          .fast_dmvnorm(densitydata, object$parameters$mean[,k], solve_t_chol_sigma= solve_t_chol_sigma_list[[k]],log=log)
      }
      mixture <- rowSums(density) ## sum(Li) 
      if ( !is.null(Sobs_activeBoundaries)) mixture <- mixture*freqs
    } else mixture <- freqs
  }
  return(mixture)
}

