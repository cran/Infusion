## fast version of dmvnorm that uses precomputed version of chol(sigma) and minimal checking
.fast_dmvnorm <- function (x, mean, solve_t_chol_sigma, log = FALSE) {
  if (is.vector(x)) dim(x) <- c(1,length(x))
  if (nrow(x)==1L) { 
    tmp <- tcrossprod(solve_t_chol_sigma, x-mean)
    rss <- sum(tmp^2)
  } else {
    tmp <- solve_t_chol_sigma %*% (t(x)-mean) # backsolve(chol_sigma, t(x) - mean, transpose = TRUE)
    rss <- colSums(tmp^2)
  }
  logdet <- attr(solve_t_chol_sigma,"logdet")
  if (is.null(logdet)) logdet <- sum(log(diag(solve_t_chol_sigma))) # is (2017) method, precomputation is not interesting 
  logretval <- logdet - 0.5 * ncol(x) * log(2 * pi) - 0.5 * rss
  names(logretval) <- rownames(x)
  if (log) { logretval } else exp(logretval)
}

.marginalize_Rmixmod <- function(jointdens, colNames, For, over) {
  nbCluster <- jointdens@nbCluster
  margdens <- jointdens
  MEAN <- margdens@parameters@mean
  colnames(MEAN) <- colNames
  margdens@parameters@mean <- MEAN[,For,drop=FALSE]
  for (clu_it in seq_len(nbCluster)) {
    COV <- margdens@parameters@variance[[clu_it]]
    colnames(COV) <- rownames(COV) <- colNames
    margdens@parameters@variance[[clu_it]] <- COV[For,For,drop=FALSE]
  }
  margdens@varNames <- For
  margdens@simuls_activeBoundaries <- NULL
  #margdens_fP@criterionValue <- margdens_fP@likelihood <- "removed for safety"  
  return(margdens)
}

.conditional_Rmixmod <- function(jointdens, fittedPars, given, expansion=Infusion.getOption("expansion")) {
  nbCluster <- jointdens@nbCluster
  conddens <- jointdens
  MEAN <- conddens@parameters@mean
  givenNames <- names(given)
  colnames(MEAN) <- colNames <- c(fittedPars, givenNames)
  For <- fittedPars 
  conddens@parameters@mean <- MEAN[,For,drop=FALSE] # resizes, but will be modified
  for (clu_it in seq_len(nbCluster)) {
    COV <- conddens@parameters@variance[[clu_it]]
    colnames(COV) <- rownames(COV) <- colNames
    sig22 <-  COV[givenNames,givenNames,drop=FALSE]
    sig12 <-  COV[For,givenNames,drop=FALSE]
    mean2 <- MEAN[clu_it,givenNames]
    conddens@parameters@mean[clu_it,] <- MEAN[clu_it,For] + sig12 %*% solve(sig22,given-mean2)
    sig11 <- COV[For,For,drop=FALSE]
    conddens@parameters@variance[[clu_it]] <- expansion* (sig11 - sig12 %*% solve(sig22,t(sig12))) 
  }
  conddens@varNames <- For
  conddens@simuls_activeBoundaries <- NULL
  #margdens_fP@criterionValue <- margdens_fP@likelihood <- "removed for safety"  
  return(conddens)
}

.conditional_mclust <- function(jointdens, fittedPars, given, expansion=Infusion.getOption("expansion")) {
  nbCluster <- jointdens$G
  conddens <- jointdens
  For <- fittedPars 
  conddens$data <- conddens$data[,For]
  conddens$parameters$variance$d <- ncol(conddens$data)
  MEAN <- conddens$parameters$mean
  givenNames <- names(given)
  colNames <- c(fittedPars, givenNames)
  conddens$parameters$mean <- MEAN[For,,drop=FALSE] # resizes, but will be modified
  sigma11 <- conddens$parameters$variance$sigma[For,For,,drop=FALSE] # resizes, but will be modified
  for (clu_it in seq_len(nbCluster)) {
    sigma_it <- conddens$parameters$variance$sigma[,,clu_it] # from single array for all clusters # drops the clu_it dimension
    sig22 <-  sigma_it[givenNames,givenNames]
    sig12 <-  sigma_it[For,givenNames]
    mean2 <- MEAN[givenNames,clu_it]
    conddens$parameters$mean[,clu_it] <- MEAN[For,clu_it] + sig12 %*% solve(sig22,given-mean2)
    sigma11[,,clu_it] <- expansion* (sigma_it[For,For] - sig12 %*% solve(sig22,t(sig12))) 
  }
  conddens$parameters$variance <- .do_call_wrap("sigma2decomp", 
                                                list(sigma=sigma11,tol=-1), # tol to avoid a simplification in the structure of $orientation that may hit a but in mclust...
                                                pack="mclust") 
  # sigma2decomp determines it is of model "VVV". We must copy this info so that the right simulation function is called by Mclust:
  conddens$modelName <- conddens$parameters$variance$modelName
  return(conddens)
}


.marginalize_mclust <- function(jointdens, colNames, For, over) {
  if (length(For)==1L) stop("mclust cannot be used on one-parameter models.") # bc sigma2decomp -> cor(1,1) fails when there is a single parameter
  nbCluster <- jointdens$G
  margdens <- jointdens
  margdens$data <- margdens$data[,For]
  margdens$parameters$variance$d <- ncol(margdens$data)
  margdens$parameters$mean <- margdens$parameters$mean[For,,drop=FALSE]
  COV <- margdens$parameters$variance$sigma # single array for all clusters 
  COV <- COV[For,For,,drop=FALSE]
  margdens$parameters$variance <- .do_call_wrap("sigma2decomp", 
                                                list(sigma=COV,tol=-1), # tol to avoid a simplification in the structure of $orientation that may hit a but in mclust...
                                                pack="mclust") 
  # Not necess useful here, but following thesamelogic as in .conditional_mclust():
  margdens$modelName <- margdens$parameters$variance$modelName
  return(margdens)
}



.get_best_clu_by_AIC <- function(cluObject) { ## cluObject: mixmodCluster object
  BICs <- unlist(lapply(cluObject@results,slot,name="criterionValue"))
  logLs <- unlist(lapply(cluObject@results,slot,name="likelihood"))
  dfs <- (2*logLs+BICs)/(log(cluObject@nbSample))
  AICs <- -2*logLs+2*dfs
  return(cluObject@results[[which.min(AICs)]])
}

.solve_t_cholfn <- function(mat) {
  resu <- solve(t(chol(mat)))
  return(structure(resu,logdet=sum(log(diag(resu)))))
}

infer_SLik_joint <- function(data, ## reference table ~ abc
                            stat.obs,
                            logLname=Infusion.getOption("logLname"), ## not immed useful
                            Simulate=attr(data,"Simulate"), ## may be NULL
                            nbCluster=Infusion.getOption("nbCluster"),
                            using=Infusion.getOption("using"),
                            verbose=list(most=interactive(), ## must be explicitly set to FALSE in knitr examples
                                         pedantic=FALSE,
                                         final=FALSE),
                            marginalize=TRUE
) {
  if ( ! is.data.frame(data)) {
    stop(paste("'object' is not a data.frame.\n Did you mean to call infer_logLs() rather than infer_Slik_joint() ?"))
  }
  if ( ! is.null( cn <- colnames(stat.obs))) {
    message("Note: 'stat.obs' should be a numeric vector, not a matrix or data.frame. Converting...")
    stat.obs <- drop(stat.obs)
    names(stat.obs) <- cn ## minimal patch so that names() can be used, not colnames()
  }
  if (!is.list(verbose)) verbose <- as.list(verbose)
  if (is.null(names(verbose))) names(verbose) <- c("most","final")[seq_len(length(verbose))]
  if (is.null(verbose$most)) verbose$most <- interactive()
  if (is.null(verbose$pedantic)) verbose$pedantic <- FALSE
  if (is.null(verbose$final)) verbose$final <- FALSE
  if (length(unique(colnames(data)))<ncol(data)) {
    stop("Some columns of the 'data' have identical names. Please correct.")
  }
  statNames <- names(stat.obs)
  for (st in statNames) {
    tabstat <- table(data[,st])
    #which_mass <- as.numeric(names(which(tabstat>1L)))
    if (any(tabstat>Infusion.getOption("nodesize"))) { # this is a projector-specific value yet we use a single value fro the slik object
      # There is <ranger object>$min.node.size but no such info in the RandomForest results.
      massvalues <- as.numeric(names(which.max(tabstat))) 
      names(massvalues) <- rep(st,length(massvalues))
      attr(stat.obs,"boundaries") <- c(attr(stat.obs,"boundaries"),massvalues)
    }
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
  if (using=="mclust") {
    if ( ! "package:mclust" %in% search()) stop("'mclust' should be loaded first.")
    if (verbose$most) cat(paste("Densities modeling:",nrow(data),"points;\n"))
    
    if (length(nbCluster$jointdens)==1L) {
      jointdens <- .do_call_wrap("densityMclust",
                                 list(data=data[,c(fittedPars,statNames)],modelNames=Infusion.getOption("mclustModel"), 
                                      G=nbCluster$jointdens,verbose=FALSE),
                                 pack="mclust")
    } else { 
      models <- vector("list",length(nbCluster$jointdens))
      for (it in seq_along(nbCluster$jointdens)) {
        models[[it]] <- .do_call_wrap("densityMclust",
                                      list(data=data[,c(fittedPars,statNames)],modelNames=Infusion.getOption("mclustModel"), 
                                           G=nbCluster$jointdens[it],verbose=FALSE),
                                      pack="mclust")
      }
      jointdens <- .get_best_mclust_by_IC(models) 
    }
    if (verbose$most) cat(paste(" joint density modeling:",jointdens$G,"clusters;\n"))
    if (verbose$pedantic) if (jointdens$G==max(nbCluster$jointdens)) message("Inferred # of clusters= max of allowed range.")
    if (marginalize) {
      pardens <- .marginalize_mclust(jointdens, colNames= c(fittedPars,statNames),
                              For=fittedPars, over=statNames)
      # defective but should be sufficient for prediction. Some garde-fou:
      #pardens@proba <- matrix(NA_real_)
    } else {
      if (length(nbCluster$pardens)==1L) {
        pardens <- .do_call_wrap("densityMclust",
                                 list(data=data[,fittedPars,drop=FALSE],modelNames=Infusion.getOption("mclustModel"), G=nbCluster$pardens,verbose=FALSE),
                                 pack="mclust")
      } else { 
        models <- vector("list",length(nbCluster$pardens))
        for (it in seq_along(nbCluster$pardens)) {
          models[[it]] <- .do_call_wrap("densityMclust",
                                        list(data=data[,fittedPars,drop=FALSE],modelNames=Infusion.getOption("mclustModel"), 
                                             G=nbCluster$pardens[it], verbose=FALSE),
                                        pack="mclust")
        }
        pardens <- .get_best_mclust_by_IC(models) 
      }
      if (verbose$most) cat(paste(" parameter modeling:",pardens$G,"clusters.\n"))
    }
    resu <- list(jointdens=jointdens,pardens=pardens)
    # pdl <- vector("list", pardens$G)
    # sigma <- resu$pardens$parameters$variance$sigma
    # for (it in seq_len(pardens$G)) pdl[[it]] <- solve(t(chol(sigma[,,it])))
    # jdl <- vector("list", jointdens$G)
    # sigma <- resu$jointdens$parameters$variance$sigma
    # for (it in seq_len(jointdens$G)) jdl[[it]] <- solve(t(chol(sigma[,,it])))
    # resu$solve_t_chol_sigma_lists <- list(pardens=pdl, jointdens=jdl) 
  } else {
    models <- .do_call_wrap("mixmodGaussianModel",list(listModels=Infusion.getOption("mixmodGaussianModel")))
    if (verbose$most) cat(paste("Densities modeling:",nrow(data),"points;\n"))
    locarglist <- list(data=data[,c(fittedPars,statNames)],nbCluster=nbCluster$jointdens, 
                       seed=Infusion.getOption("mixmodSeed") , models=models)
    checkDegenerate <- cov(locarglist$data)
    eigvals <- eigen(checkDegenerate, only.values = TRUE)$values
    if (any(abs(eigvals<1e-14))) warning(paste(
      c("The covariance matrix of the (parameters,statistics) table seems singular,\n",
        "implying linear relationships between the variables. Problems will likely happen.\n",
        "Redundant variables should be eliminated."
        )),immediate. = TRUE)
    if (FALSE) { ## older code; still useful for devel
      jointdens <- try(.do_call_wrap("mixmodCluster",locarglist),silent = TRUE)
      jointdens <- .get_best_clu_by_AIC(jointdens)
    }
    jointdens <- try(.do_call_wrap(".densityMixmod",c(locarglist,list(stat.obs=stat.obs))),silent = TRUE) 
    if (inherits(jointdens,"try-error")) {
      nbCluster <- eval(Infusion.getOption("nbCluster"))
      if ( ! identical(nbCluster, locarglist$nbCluster)) {
        locarglist$nbCluster <- nbCluster
        jointdens <- .do_call_wrap(".densityMixmod",c(locarglist,list(stat.obs=stat.obs)))
      } else stop(jointdens)
    } ## else .densityMixmod has already run .get_best_clu_by_AIC()  
    # plotCluster(jointdens,data=locarglist$data,variable1="theta_p",variable2="theta") # var1: stat (prediction of projection); var2: actual param
    if (verbose$most) cat(paste(" joint density modeling:",jointdens@nbCluster,"clusters;\n"))
    if (verbose$pedantic) if (jointdens@nbCluster==max(nbCluster$jointdens)) message("Inferred # of clusters= max of allowed range.")
    if (marginalize) {
      pardens <- .marginalize_Rmixmod(jointdens, colNames= c(fittedPars,statNames),
                              For=fittedPars, over=statNames)
      # defective but should be sufficient for prediction. Some garde-fou:
      pardens@proba <- matrix(NA_real_)
      pardens@parameters@nbFreeParam <- NA_integer_
    } else {
      locarglist <- list(data=data[,fittedPars,drop=FALSE], nbCluster=nbCluster$pardens, 
                         seed=123 , models=models)
      pardens <- try(.do_call_wrap("mixmodCluster",locarglist),silent = TRUE)
      if (inherits(pardens,"try-error")) {
        locarglist$nbCluster <- eval(Infusion.getOption("nbCluster"))
        pardens <- do.call("mixmodCluster",locarglist)
      }                               
      if (length(nbCluster)==1L) {
        pardens <- pardens@results[[1L]]
      } else { pardens <- .get_best_clu_by_AIC(pardens) }
      if (verbose$most) cat(paste(" parameter modeling:",pardens@nbCluster,"clusters.\n"))
    }
    #plotCluster(pardens,data=data[,fittedPars]) to plot @results[[1L]] which is a 'MixModResults', not a 'mixmodCluster' object.
    # Rmixmod::plot(<mixmodCluster object>) with Rmixmod:: to plot from any envir, not only the global one
    resu <- list(jointdens=jointdens,pardens=pardens)
    resu$solve_t_chol_sigma_lists <- list(
      pardens=lapply(resu$pardens@parameters["variance"], .solve_t_cholfn),
      jointdens=lapply(resu$jointdens@parameters["variance"], .solve_t_cholfn)
    ) 
  }
  attr(resu,"EDFstat") <- "[see this string in infer_SLik_joint()]" ## 
  resu$logLs <- structure(data,stat.obs=stat.obs,Simulate=Simulate) ## as in infer_surface.logLs
  # attr(resu,"callfn") <- "infer_SLik_joint"
  resu$projectors <- attr(data,"projectors")
  resu$`Infusion.version` <- packageVersion("Infusion")
  resu$colTypes <- list(allPars=allPars, ## keeps the order of the columns
                                fittedPars=fittedPars,
                                fixedPars=fixedPars,
                                logLname=logLname,
                                statNames=statNames)
  resu$lower <- apply(data[,fittedPars,drop=FALSE],2L,min) ## really this
  resu$upper <- apply(data[,fittedPars,drop=FALSE],2L,max) # ... used by MSL -> optim
  if (any(is.na(c(resu$lower,resu$upper)))) stop("NA(s) in c(lower,upper))") 
  resu$LOWER <- attr(data,"LOWER") # used in many places
  if (is.null(resu$LOWER)) resu$LOWER <- resu$lower
  resu$UPPER <- attr(data,"UPPER") # ... used by confintAll
  if (is.null(resu$UPPER)) resu$UPPER <- resu$upper
  attr(resu,"Infusion.version") <- packageVersion("Infusion")
  class(resu) <- c("SLik_j",class(resu))
  predEDF <- predict(resu,newdata=data[,fittedPars,drop=FALSE])
  wm <- which.max(predEDF)
  resu$optrEDF <- list(par=data[wm,fittedPars,drop=FALSE],value=predEDF[wm])
  resu$raw_data <- attr(data,"raw_data")
  return(resu)
}


## iwas nfer_SLik_joint -> predict -> predict.SLik_j -> predict.MixmodResults
# except that now                                    -> .pointpredict.Rmixmod -> predict.dMixmod
predict.MixmodResults <- function(object, newdata,log=TRUE, solve_t_chol_sigma_list,...) {
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
      if (is.vector(normalizedlogs)) dim(normalizedlogs) <- c(1L, length(normalizedlogs))
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

.pointpredict.Rmixmod <- function(point, object, tstat.obs, log, which, ...)  {
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

.predict_SLik_j_mclust <- function(object, newdata, tstat.obs, log, which)  {
  if (which!="parvaldens") jointvaldens <- predict(object$jointdens, newdata=cbind(data.frame(newdata),tstat.obs)) # ! order of elements in newdata must be that of fittedPars as in object$jointdens$data
  if (which=="jointvaldens") return(jointvaldens)
  parvaldens <- predict(object$pardens,newdata=newdata) 
  if (which=="parvaldens") return(parvaldens)
  if (log) {
    jointvaldens <- pmax(jointvaldens, .Machine$double.xmin) # F I X M E for the value
    condvaldens <- log(jointvaldens) - log(parvaldens)
  } else condvaldens <- jointvaldens/parvaldens
  if (any(is.infinite(condvaldens))) warning("any(is.infinite(condvaldens)) is TRUE: expect an error.") 
  return(condvaldens)
}

predict.SLik_j <- function (object, 
                            newdata, ## requests new fittedPars values! 
                            # newdata=object$logLs[,object$colTypes$fittedPars,drop=FALSE], 
                            log=TRUE, 
                            which="condvaldens",
                            tstat= t(attr(object$logLs,"stat.obs")),
                            ...) {
  if (is.null(nrow(newdata)) ) newdata <- t(as.matrix(newdata)) # as.matrix keeps names
  if (is.null(colnames(newdata))) colnames(newdata) <- object$colTypes$fittedPars
  # not useful bc predict.MixmodResults does not handle it:
  #newdata <- cbind(newdata, t(replicate(nrow(newdata),attr(object$logLs,"stat.obs"))))
  if (inherits(object$jointdens,"densityMclust")) {
    .predict_SLik_j_mclust(object=object, newdata=newdata, tstat.obs=tstat, log=log, which=which)
  } else {
    apply(newdata,1L, .pointpredict.Rmixmod, object=object, tstat.obs=tstat, log=log, which=which)
  }
}

confint.SLik_j <- function(object, parm, ## parm is the parameter which CI is sought 
                         level=0.95, verbose=interactive(),fixed=NULL,which=c(TRUE,TRUE),...) {
  .confintAll(object=object, parm=parm, ## parm is the parameter which CI is sought 
             givenmax = object$MSL$maxlogL,
             level= - qchisq(level,df=1)/2, ## df=1 for 1D profile; /2 to compare to logLR rather than 2 logLR  
             verbose=verbose,fixed=fixed,which=which,...)
}

refine.SLik_j <- function(object,method=NULL,...) {
  if (is.null(method)) method <- "mixmodCluster" ## no clear effect but distinct from SLik case
  refine.default(object, surfaceData=object$logLs, method=method, ...)
}

print.SLik_j <-function(x,...) {summary.SLik_j(x,...)} 

.rparam_from_SLik_j <- function(object, 
                                frac, 
                                fittedPars=NULL,level) {
  if (is.null(fittedPars)) fittedPars <- object$colTypes$fittedPars
  CI_LRstat <- qchisq(level,df=length(fittedPars)) 
  size_first_iter <- attr(object$logLs,"n_first_iter")
  if (is.null(size_first_iter)) size_first_iter <- length(which(object$logLs$cumul_iter==1)) ## FR->FR temporary back compat
  prev_n_iter <- max(object$logLs$cumul_iter)
  ceil_size <- target_size <- frac*max(size_first_iter/2, size_first_iter*(prev_n_iter+1L)-nrow(object$logLs))
  if ( ! is.null(freq_good <- attr(object$logLs,"freq_good")$uniform)) ceil_size <- ceiling(ceil_size/freq_good)
  if (inherits(object$jointdens,"densityMclust")) {
    trypoints <- do.call("sim", c(object$jointdens[c("modelName", "parameters")],list(n=ceil_size)))
    trypoints <- trypoints[,1L+seq_len(length(fittedPars)),drop=FALSE]  
  } else {
    trypoints <- .simulate.MixmodResults(object$jointdens, nsim=1L, size=ceil_size, drop=TRUE) ## simulates joint
    trypoints <- trypoints[,seq_len(length(fittedPars)),drop=FALSE] # sample in marginal distrib of params from sample in jointdens
  } ## : samples the joint and keep only the params => distrib of trypoints=current distr of param (pardens)
  colnames(trypoints) <- fittedPars
  prior <- predict(object,trypoints,which="parvaldens") ## slow
  logLs <- predict(object,trypoints,which="condvaldens") ## slow
  criterion <- logLs-prior ## to *cancel* the fact that we sampled according to pardens
  threshold <- max(criterion)-CI_LRstat ## defines a slice of criterion values
  regul <- pmax(threshold-CI_LRstat, ## sets a minimal weight=exp(-CI_LRstat)
                pmin(threshold,criterion)) ## sets a maximal weight = exp(0)
  weights <- exp(regul-threshold) 
  if (FALSE) { ## ring sampling. yields poor clustering.
    locsd <- 1/5
    weights <- weights*sapply(criterion,function(v) max(dnorm(v-c(0,CI_LRstat/2),sd = locsd)))/(0.3989422804014327*locsd)
  }
  good <- which(runif(n=length(weights))<weights) ## samples uniformly in the current fitted distrib of params with logL>threshold
  freq_good <- length(good)/ceil_size
  if (freq_good < 1/5) {
    good <- order(weights,decreasing=TRUE)[seq(ceil_size/5)]
    freq_good <- 1/5 ## to avoid explosion of # of sampled points
  } else if (length(good)>target_size) {
    good <- good[sample(length(good),size = target_size)]
  } # may stay a few points below target, rather than add points that were not retained bc of low weights (it's better to wait next iteration)
  trypoints <- trypoints[good,,drop=FALSE]
  trypoints <- (cbind(trypoints,object$colTypes$fixedPars)) ## add fixedPars for simulation
  trypoints <- trypoints[,object$colTypes$allPars,drop=FALSE] ## column reordering
  trypoints <- data.frame(trypoints)
  for (vv in fittedPars) trypoints <- trypoints[trypoints[,vv]>object$LOWER[vv],,drop=FALSE]
  for (vv in fittedPars) trypoints <- trypoints[trypoints[,vv]<object$UPPER[vv],,drop=FALSE]
  return(list(trypoints=trypoints, freq_good=freq_good))
}


.rparam_from_SLik_j_post <- function(object , frac) {
  if (frac==0) return(list(trypoints=NULL, freq_good=NA))
  size_first_iter <- attr(object$logLs,"n_first_iter")
  if (is.null(size_first_iter)) size_first_iter <- length(which(object$logLs$cumul_iter==1)) ## FR->FR temporary back compat
  prev_n_iter <- max(object$logLs$cumul_iter)
  ceil_size <- target_size <- frac*max(size_first_iter/2, size_first_iter*(prev_n_iter+1L)-nrow(object$logLs))
  if ( ! is.null(freq_good <- attr(object$logLs,"freq_good")$posterior)) ceil_size <- ceiling(ceil_size/freq_good)
  fittedPars <- object$colTypes$fittedPars
  if (inherits(object$postdens,"densityMclust")) {
    trypoints <- do.call("sim", c(object$postdens[c("modelName", "parameters")],list(n=ceil_size)))
    trypoints <- trypoints[,1L+seq_len(length(fittedPars)),drop=FALSE]
  } else {
    trypoints <- .simulate.MixmodResults(object$postdens, nsim=1L, size=ceil_size, drop=TRUE)
  } ## : samples the joint and keep only the params => distrib of trypoints=current distr of param (pardens)
  colnames(trypoints) <- fittedPars
  trypoints <- (cbind(trypoints,object$colTypes$fixedPars)) ## add fixedPars for simulation
  trypoints <- trypoints[,object$colTypes$allPars,drop=FALSE] ## column reordering
  trypoints <- data.frame(trypoints)
  for (vv in fittedPars) trypoints <- trypoints[trypoints[,vv]>object$LOWER[vv],,drop=FALSE]
  for (vv in fittedPars) trypoints <- trypoints[trypoints[,vv]<object$UPPER[vv],,drop=FALSE]
  freq_good <- nrow(trypoints)/ceil_size
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

# this does not call the process-simulating function. Instead, it performs a sort of parametric bootstrap from the Gaussian mixture model.
.boot.SLik_j <- function(object,nsim=2L, force=FALSE, verbose = TRUE, seed=NULL, 
                         parent_cores_info=NULL, # defined cluster by parent *.boot.SLik_j* calls
                         nb_cores=NULL, packages=attr(object$logLs,"packages"), env=attr(object$logLs,"env") # not required if parent_cores_info is provided.
                         ) {
  if (!force && nsim<2L) stop("'nsim' must be > 1")
    if (inherits(object$jointdens,"densityMclust")) {
    bootrepls <- replicate(nsim, do.call("sim", c(object$jointdens[c("modelName", "parameters")],
                                                  list(n=nrow(object$logLs))))[,-1L],simplify=FALSE)
    nbCluster <- list(jointdens=object$jointdens$G,
                      pardens=object$pardens$G)
  } else {
    bootrepls <- .simulate.MixmodResults(object$jointdens, nsim=nsim, size=nrow(object$logLs),
                                         drop=FALSE) 
    nbCluster <- list(jointdens=object$jointdens@nbCluster,
                      pardens=object$pardens@nbCluster) # single values so that this exposes to clustering failure.
  } ##  marginal=current distr of param
  
  boo <- object$logLs[,with(object$colTypes,c(fittedPars,statNames))]
  it <- 0L
  prevmsglength <- 0L
  resu <- vector("list")
  boo_SLik_joint <- function(simul) {
    it <<- it+1
    locmess <- paste("boostrap replicate",it,"of",nsim,"   ")
    if (verbose) {prevmsglength <<- .overcat(locmess,prevmsglength = prevmsglength)}
    boo[] <- simul
    densv <- suppressWarnings( ## suppress warnings for clustering failure
      try(infer_SLik_joint(boo,stat.obs=attr(object$logLs,"stat.obs"), nbCluster=nbCluster,
                           verbose=list(most=FALSE,final=FALSE)),silent=TRUE))
    if (inherits(densv,"try-error")) {
      return(NULL)
    } else return(densv)
  }
  if (is.null(parent_cores_info)) {
    cores_info <- .init_cores(nb_cores=nb_cores)
    cl <- cores_info$cl
    if ( ! is.null(cl)) {
      parallel::clusterExport(cl, list("packages"),envir=environment()) ## passes the list of packages to load
      # Infusion not leaded in child process !
      #parallel::clusterExport(cl, list("nRealizations"),envir=environment()) 
      #parallel::clusterCall(cl,fun=Infusion.options,nRealizations=nRealizations)
      if ( ! is.null(env)) parallel::clusterExport(cl=cl, as.list(ls(env)),envir=env)
    }
  } else cores_info <- parent_cores_info
  resu <- pblapply(bootrepls,boo_SLik_joint,cl = cores_info$cl)
  #resu <- lapply(bootrepls,boo_SLik_joint)
  whichvalid <- which( ! sapply(resu,is.null))
  resu <- resu[whichvalid]
  while (length(resu)< nsim) {
    message(paste(nsim-length(resu),"failed bootstrap replicate(s); drawing sample(s) again..."))
    moreresu <- .boot.SLik_j(object, nsim=nsim-length(resu), force=TRUE, verbose=verbose, parent_cores_info=cores_info)
    resu <- c(resu,moreresu)
  }
  if (is.null(parent_cores_info)) .close_cores(cores_info)
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
# .plotClusters <- function(object,which,...) {
#   plotCluster(object[[which]],data=object$logLs,...)
#}

plot.MixmodResults <- function(x,...) {message('Use plotCluster(.,data=<slik object>$logLs,variable1=.,variable2=.)\n to plot a "MixModResults" object')}
